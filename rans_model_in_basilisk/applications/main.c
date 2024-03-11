#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wmultichar"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
/**
# Breaking Stokes wave

We solve the two-phase Navier--Stokes equations using a
momentum-conserving transport of each phase. Gravity is taken into
account using the "reduced gravity approach". */
//#include <unistd.h>
//#include <stdlib.h>
//#include <stdio.h>
//#include <error.h>
//#include <argp.h>

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "reduced.h"
#include "tag.h"

#include "sander/output_htg.h"
#include "oystein/adapt_wavelet_leave_interface.h"
#include "oystein/output_surfaces.h"

int kVerbose = 0;
int kAbort = 0;
#include "json_utils.h"
#include "params.h"
#include "stokes_t.h"
#include "focus_t.h"

#include "rsm_common.h"

/************************************************************
Global variables
************************************************************/
extern struct WaveProperties kWaveData;
extern struct MultidirectionInletProfile kMultiwaveInlet;
extern struct WaveGaugeProfile kWaveGauge;
struct rsm_JSONFile_t kInputJson = {NULL,0,NULL};
struct rsm_JSONFile_t kMultiWaveInletJson = {NULL,0,NULL};
struct rsm_JSONFile_t kWaveGaugeJson = {NULL,0,NULL};

#ifdef RSM_STRESS_OMEGA
  extern struct StressOmegaConstants kTurbConstants;
#endif
#ifdef RSM_KOMEGA
  extern struct KOmegaConstants kTurbConstants;
#endif

scalar wall_mask[]; // indicator field for refinement
vector gp[];
scalar dp[];
/************************************************************
Static variables and forward declarations
************************************************************/
#define M_PI (3.1415926)
// NOTE: These static variables are determined during compilation hence cannot be modified in runtime.
static const double SURFACE_OUTPUT_PERIOD = 0.01; // second
static const double OUTPUT_PERIOD=0.1; // second
static const double AUTOSAVE_PERIOD=1.0; // second
static const double CHECKPOINT_PERIOD=10.0; // second
static const double MAXIMUM_RUNTIME = 300.0; // second
static const double SUBGRID_DROPLET_REMOVAL_PERIOD=1.0; // second

void set_domain_size();
void set_kinematics();
void set_wall();
double ramp_frac();
double set_alpha(double x, double y, double t, double del);
/************************************************************
Boundary conditions
************************************************************/
void set_essential_boundary_conditions(){
  // Left boundary (wave inflow)
  #ifndef NO_INLET
  u.n[left]  = dirichlet(u_x(x,y,t+kWaveData.t0) * ramp_frac()) * (f[] > 0.1 ? 1.0:0.0);
  u.t[left]  = dirichlet(u_y(x,y,t+kWaveData.t0) * ramp_frac()) * (f[] > 0.1 ? 1.0:0.0);
  f[left]    = set_alpha(x,y,t+kWaveData.t0,Delta);
  #endif

  // Top boundary
  u.n[top] = neumann(0);
  p[top] = neumann(0); // this is a neat little trick to avoid escalating back-circulating flows in the top boundary.

  // bottom boundary
  u.n[bottom] = dirichlet(0);
  u.t[bottom] = dirichlet(0);

  // Additional refinement indicator
  wall_mask[bottom] = dirichlet(1);
}

#ifdef RSM_STRESS_OMEGA
void set_rsm_boundary_conditions(){
  // Left boundary (wave inflow)
#ifndef NO_INLET
  rhoe[left] = dirichlet(kTurbConstants.omega_0);
  Rxx[left] = dirichlet (-2.0/3.0*kTurbConstants.k_0);
  Rxy[left] = dirichlet (0);
  Rxz[left] = dirichlet (0);
  Ryy[left] = dirichlet (-2.0/3.0*kTurbConstants.k_0);
  Ryz[left] = dirichlet (0);
  Rzz[left] = dirichlet (-2.0/3.0*kTurbConstants.k_0);
#endif

  // Top boundary
  rhoe[top] = dirichlet(kTurbConstants.omega_0);
  Rxx[top] = dirichlet (-2.0/3.0*kTurbConstants.k_0);
  Rxy[top] = dirichlet (0);
  Rxz[top] = dirichlet (0);
  Ryy[top] = dirichlet (-2.0/3.0*kTurbConstants.k_0);
  Ryz[top] = dirichlet (0);
  Rzz[top] = dirichlet (-2.0/3.0*kTurbConstants.k_0);

  // bottom boundary
  rhoe[bottom] = dirichlet ( 10.0*(6.0*kWaveData.mu1/kWaveData.rho1) / (0.075*sq((kWaveData.size[1]-kWaveData.size[0])/(1<<kWaveData.adapt_refine_max))) );
  Rxx[bottom] = dirichlet (0);
  Rxy[bottom] = dirichlet (0);
  Rxz[bottom] = dirichlet (0);
  Ryy[bottom] = dirichlet (0);
  Ryz[bottom] = dirichlet (0);
  Rzz[bottom] = dirichlet (0);
}
#else
#ifdef RSM_KOMEGA
void set_rsm_boundary_conditions(){
  // Left boundary (wave inflow)
#ifndef NO_INLET
  rhoe[left] = dirichlet(kTurbConstants.omega_0);
  rhok[left] = dirichlet (kTurbConstants.k_0);
#endif

  // Top boundary
  rhoe[top] = dirichlet(kTurbConstants.omega_0);
  rhok[top] = dirichlet (kTurbConstants.k_0);

  // bottom boundary
  rhoe[bottom] = dirichlet ( 10.0*(6.0*kWaveData.mu1/kWaveData.rho1) / (0.075*sq((kWaveData.size[1]-kWaveData.size[0])/(1<<kWaveData.adapt_refine_max))) );
  rhok[bottom] = dirichlet (kTurbConstants.kMin_);
}
#else
void set_rsm_boundary_conditions(){}
#endif
#endif

/************************************************************
Argument Parser
************************************************************/
const char *argp_program_version =
  "waveflume v1.0";
const char *argp_program_bug_address =
  "<yuxuan.liu@eng.ox.ac.uk>";

/* Program documentation. */
static char doc[] =
  "Waveflume -- a cfd program to simulate ocean wave dynamics;\n\
options\
\vThis program is implemented under Basilisk framework;\n\
 url to Basilisk official website.";

/* A description of the arguments we accept. */
static char args_doc[] = "";

/* Keys for options without short-options. */
#define OPT_ABORT  1            /* abort */
/*
// The options we understand. //
static struct argp_option options[] = {
  // System configuration
  {"abort",  OPT_ABORT, 0,       0, "Abort the program early" },
  {"verbose",  'v', "INT",       0, "Produce verbose output, level 0 - 3" },
  {"config", 'i', "FILE", 0, "Path to configuration file"},
  {"multi_inlet_config", 'm', "FILE", 0, "Path to multiwave inlet configuration file"},

  { 0 }
};

// Parse a single option. //
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  // Get the input argument from argp_parse, which we
  // know is a pointer to our arguments structure. //
  //struct WaveProperties *arguments = ((void**)state->input)[0];

  switch (key)
  {
    // System configuration
    case OPT_ABORT:
      kAbort = 1;
      break;

    case 'v':
	    kVerbose = arg?atoi(arg):1;
	    break;

    case 'i':
      kInputJson.path = arg;
      break;

    case 'm':
      kMultiWaveInletJson.path = arg;
      break;

    default:
      return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

// Our argp parser. //
static struct argp argp = { options, parse_opt, args_doc, doc };
*/
void parse_global_vars_by_json(struct rsm_JSONFile_t *jf)
{
    read_int(&kVerbose, "verbose", jf->json);
}

/************************************************************
The program takes optional arguments which are the level of
refinement, steepness and Reynolds numbers.
************************************************************/
int main (int argc, char * argv[])
{
  /* Parse our arguments; every option seen by parse_opt will be
	 reflected in arguments. */
  //argp_parse (&argp, argc, argv, 0, 0, NULL);
  if(argc > 1)
    kInputJson.path = argv[1];

  if(argc > 2)
    kMultiWaveInletJson.path = argv[2];

  if(argc > 3)
    kWaveGaugeJson.path = argv[3];
  
  if(argc > 4)
    kAbort = atoi(argv[4]);

  // Read json input file
  if(!kInputJson.path)
    kInputJson.path = "./config.json"; // set to a default path if not initialized by other subroutines...
  update_configuration_file(&kInputJson);
  parse_global_vars_by_json(&kInputJson);

  // update multidirectional wave inlet configuration
  if(!kMultiWaveInletJson.path)
    kMultiWaveInletJson.path = "./wave_profile.json"; // set to a default path if not initialized by other subroutines...
  update_configuration_file(&kMultiWaveInletJson);

  // update wave gauge configuration
  if(!kWaveGaugeJson.path)
    kWaveGaugeJson.path = "./wave_gauge.json"; // set to a default path if not initialized by other subroutines...
  update_configuration_file(&kWaveGaugeJson);

  refresh_waveinput_by_json(&kInputJson, &kMultiWaveInletJson);
  refresh_multiwave_inlet_by_json(&kMultiWaveInletJson);
  refresh_wave_gauge_by_json(&kWaveGaugeJson);

  // update turbulence model coeffs
  refresh_rsm_model_by_json(&kInputJson);

#if _MPI  
  MPI_Barrier(MPI_COMM_WORLD);  
#endif

  if (pid() == 0 && kAbort == 0)
	  system ("mkdir -p postProcessing");

  /**
  The domain is a cubic box centered on the origin and of length
  $L0=1$, periodic in the x-direction. */
  set_domain_size();
#ifdef NO_INLET
  periodic (right);
#endif
#if dimension == 3
  periodic (back);
#endif

  /**
  Here we set the densities and viscosities corresponding to the
  parameters above. */
  rho1 = kWaveData.rho1;
  rho2 = kWaveData.rho2;
  mu1 = kWaveData.mu1; //using wavelength as length scale
  mu2 = kWaveData.mu2;
  G.y = -kWaveData.g;


  /**
  When we use adaptive refinement, we start with a coarse mesh which
  will be refined as required when initialising the wave. */
  
#if TREE  
  N = 1 << kWaveData.refine_init;
#else
  N = 1 << kWaveData.adapt_refine_max;
#endif
  DT = kWaveData.dtmax;
  CFL = kWaveData.cfl;

  if (pid()==0)
    PrintWaveProperties(&kWaveData);

  run();
}

/**
## Initial conditions

We either restart (if a "restart" file exists), or initialise the wave
using the third-order Stokes wave solution. */
event defaults ( i = 0 ) {
  set_essential_boundary_conditions();
  set_rsm_boundary_conditions();

  // Set the previous timestep to DT
  face vector uf_dummy[];
  foreach_face()
    uf_dummy.x[] = 1e-30;
  double dt_dummy = timestep (uf_dummy, 1000);
}

event init (i = 0) {
  fprintf(stdout,"###Main initialization subroutine ###\n");

  if (!restore ("restart")) {
    int count = kWaveData.adapt_refine_max - kWaveData.refine_init + 1;
    do {	
	set_kinematics();
	
	// modify indicator field
	set_wall();
	
	if(count-- < 0){
	  break;
	}
  fprintf(stdout,"Init iteration %d\n", count);

    }
    /**
    On trees, we repeat this initialisation until mesh adaptation does
    not refine the mesh anymore. */
#if TREE  
    while (
      #if dimension == 2
      adapt_wavelet_leave_interface((scalar *){u, wall_mask},{f},
      (double[]){kWaveData.adapt_utol,kWaveData.adapt_utol,0.001}, 
      kWaveData.adapt_refine_max, kWaveData.adapt_refine_max, kWaveData.adapt_refine_min,1).nf
      #elif dimension == 3
      adapt_wavelet_leave_interface((scalar *){u, wall_mask},{f},
      (double[]){kWaveData.adapt_utol,kWaveData.adapt_utol,kWaveData.adapt_utol,0.001}, 
      kWaveData.adapt_refine_max, kWaveData.adapt_refine_max, kWaveData.adapt_refine_min,1).nf
      #endif
    );
#else
    while (0);
#endif
  }
}

#ifndef NO_CONFIG_REFRESH
event watch_configuration_file(i++) {
  if(update_configuration_file_when_modified(&kInputJson) == 0x1)
  {
    rsm_print_jsonfile(&kInputJson);
    parse_global_vars_by_json(&kInputJson);
    refresh_waveinput_by_json(&kInputJson, &kMultiWaveInletJson);
    refresh_rsm_model_by_json(&kInputJson);
  }
}
#endif

event remove_bubbles0(t+=SUBGRID_DROPLET_REMOVAL_PERIOD, last) {
  remove_droplets(f);
  fprintf (stderr, "[%d] %g %d remove droplets\n", pid(), t, i);
}


event stability (i++) {
  CFL = kWaveData.cfl;
}


event logfile (i+=100;t <= MAXIMUM_RUNTIME) {
	int world_rank = -1; // the rank of the process
# if _MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#endif
	stats s = statsf (f);
	fprintf (stderr, "[%d] %g %d %g %g %.8f\n", world_rank, t, i, s.min, s.max, s.sum);
}

/**
Save a checkpoint. */
event autosave(t += AUTOSAVE_PERIOD) {
foreach(){
  dp[] = Delta;
  gp.x[] = x;
  gp.y[] = y;
  #if dimension > 2
    gp.z[] = z;
  #endif
}
#if _MPI  
  MPI_Barrier(MPI_COMM_WORLD);  
#endif
  dump(file="restart", unbuffered=false);
#if _MPI  
  MPI_Barrier(MPI_COMM_WORLD);  
#endif
  fprintf (stderr, "%d : Dump restart file.\n", i);
  if (kAbort){
#if _MPI  
    MPI_Finalize();
    exit(0);
#endif
  }
}

event checkpoint(t += CHECKPOINT_PERIOD) {
#if _MPI  
  MPI_Barrier(MPI_COMM_WORLD);  
#endif
  char fname[50];
  sprintf(fname, "restart.%f.ckpt", t);
  dump(file=fname, unbuffered=false);
#if _MPI  
  MPI_Barrier(MPI_COMM_WORLD);  
#endif
  fprintf (stderr, "%d : Dump checkpoint file.\n", i);
}

/**
Write to hypretree format. */
event write2htg (t += OUTPUT_PERIOD) {
  scalar vort[], p_rgh[];
  vorticity (u, vort);

  vector U[]; // to slove the compatibility issue of restore/dump
  foreach()
    foreach_dimension()
      U.x[] = u.x[];
  boundary((scalar *){U});

  char fname[50];
  sprintf(fname, "result.%06d", i);
  scalar* output_scalars = {f, p, vort};
  rsm_model_output(&output_scalars);
  output_htg(output_scalars,(vector *){U}, ".", fname, i, t);

  fprintf (stderr, "write output to %s\n", fname);
}

/**
Write ply and csv output. */
event outputInterface(t += SURFACE_OUTPUT_PERIOD) {
  char resultname[100];
  sprintf( resultname, "postProcessing/results_%4.2f.dat", t );
  FILE * fp = fopen(resultname, "w");

  scalar* output_scalars = {p};
  rsm_model_output_compact(&output_scalars);
  reconstruct_interface(f,kWaveGauge.rays, kWaveGauge.num_probes, fp,output_scalars,(vector *){u});

  fclose (fp);
}

event finishSimulation(i++;t<=MAXIMUM_RUNTIME) {
#if _MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if ( t >= (kWaveData.tmax - kWaveData.t0) ){
#if _MPI  
    MPI_Finalize();
#endif
    exit(0);
  }
}

//------------------DIAGNOSTICS---------------------//
//scalar dissrate[];

double sech(double qval) {
  return 1.0/cosh(qval);
}

double maxv (scalar a) {
  double maxi = - 1.0e100;

//Calculate the maximum quantity in the water.

  foreach (reduction(max:maxi)){
    if (fabs(a[]*f[]) > maxi)
      maxi = fabs(a[]*f[]);
  }
  return maxi;
}

/*Define functions for determining kinetic and potential energy*/
int dissipation_rate (vector u, double* rates)
{
  double rateWater = 0.0;
  double rateAir = 0.0;
  foreach (reduction (+:rateWater) reduction (+:rateAir)) {
    double dudx = (u.x[1] - u.x[-1])/(2.0*Delta);
    double dudy = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]   - u.y[-1]  )/(2.*Delta);
    double dvdy = (u.y[0,1] - u.y[0,-1])/(2.0*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.0*Delta);
    double SDeformxx = dudx;
    double SDeformxy = 0.5*(dudy + dvdx);
    double SDeformxz = 0.5*(dudz + dwdx);
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy;
    double SDeformyz = 0.5*(dvdz + dwdy);
    double SDeformzx = SDeformxz;
    double SDeformzy = SDeformyz;
    double SDeformzz = dwdz; 
    double sqterm = 2.0*dv()*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
			      sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
			      sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz));
    double localRate = kWaveData.mu1/rho[]*f[]*sqterm;
    //drate[] = localRate;
    rateWater += localRate; //water
    rateAir += kWaveData.mu2/rho[]*(1.0-f[])*sqterm; //air
  }
  //fprintf (fix, "\n");
  //fprintf (fiy, "\n");
  rates[0] = rateWater;
  rates[1] = rateAir;
  return 0;
}

event graphs (t+=SURFACE_OUTPUT_PERIOD, last) {
  scalar umag[];
  static FILE * fpwater = fopen("budgetWater.dat", "w");
  static FILE * fpair = fopen("budgetAir.dat", "w");
  double ke = 0., gpe = 0.;
  double keAir = 0., gpeAir = 0.;
  foreach(reduction(+:ke) reduction(+:gpe) 
	  reduction(+:keAir) reduction(+:gpeAir)) {
    double norm2 = 0.;
    foreach_dimension()
      norm2 += sq(u.x[]);
    ke += rho[]*norm2*f[]*dv()*(1.0-0);
    keAir += rho[]*norm2*(1.0-f[])*dv()*(1.0-0);
    gpe += rho[]*kWaveData.g*y*f[]*dv()*(1.0-0);
    gpeAir += rho[]*kWaveData.g*y*(1.0-f[])*dv()*(1.0-0);
    umag[] = sqrt(sq(u.x[]) + sq(u.y[]));
  }
  double rates[2];
  dissipation_rate(u, rates);

  double dissWater = rates[0] * rho1;
  double dissAir   = rates[1] * rho2;

  if (i == 0) {
    fprintf (fpwater, "t ke gpe dissipation\n");
    fprintf (fpair, "t ke gpe dissipation\n");
  }
  double maxs = maxv(umag);
  double gpe0    =-((kWaveData.size[1]-kWaveData.size[0])*kWaveData.rho1*kWaveData.g*0.5*sq(kWaveData.waterlevel - kWaveData.size[2]));
  double gpeAir0 = ((kWaveData.size[1]-kWaveData.size[0])*kWaveData.rho2*kWaveData.g*0.5*sq(kWaveData.size[3] - kWaveData.waterlevel));

  fprintf (fpwater, "%e %e %e %e\n",
	   t, ke/2., gpe - gpe0, dissWater);
  fprintf (fpair, "%e %e %e %e\n",
	   t, keAir/2., gpeAir - gpeAir0, dissAir);
  fprintf (ferr, "%e %e %e %e %e %e\n", dt, 
	   t, ke/2. + gpe - gpe0, ke/2., gpe - gpe0, dissWater);
}

/**
## Mesh adaptation

On trees, we adapt the mesh according to the error on volume fraction
and velocity. */
#if TREE
event adapt (i++) {
  #if dimension == 2
	adapt_wavelet_leave_interface((scalar *){u, wall_mask},{f},
  (double[]){kWaveData.adapt_utol,kWaveData.adapt_utol,0.001}, 
  kWaveData.adapt_refine_max>13?13:9, kWaveData.adapt_refine_max, kWaveData.adapt_refine_min, 1);
  #elif dimension == 3	
  adapt_wavelet_leave_interface((scalar *){u, wall_mask},{f},
  (double[]){kWaveData.adapt_utol,kWaveData.adapt_utol,kWaveData.adapt_utol,0.001}, 
  kWaveData.adapt_refine_max>9?9:7, kWaveData.adapt_refine_max, kWaveData.adapt_refine_min, 1);
  #endif
}
#endif

/* Utility functions */
void set_domain_size(){

	double lx = kWaveData.size[1]-kWaveData.size[0];
	double lz = kWaveData.size[3]-kWaveData.size[2];

	init_grid (1 << kWaveData.refine_init);

#if dimension == 2
	if (lx >= lz){
		size (lx);
		origin (kWaveData.size[0], kWaveData.size[2]); // move origin
	}
	else {
		size (lz);
		origin (kWaveData.size[0], kWaveData.size[2]); // move origin
	}
#elif dimension == 3
	if (lx >= lz){
		size (lx);
		origin (kWaveData.size[0], kWaveData.size[2],-lx/2); // move origin
	}
	else {
		size (lz);
		origin (kWaveData.size[0], kWaveData.size[2],-lz/2); // move origin
	}
#endif
}

void set_kinematics(){
    fraction (f, -y + wave(x,y,t+kWaveData.t0) * ramp_frac());

    foreach(){
      double eta = wave(x,y,t+kWaveData.t0) +  (kWaveData.size[1]-kWaveData.size[0])/(1<<kWaveData.adapt_refine_max);
      double d0 = x;
      double d1 = y > eta ? eta : y;
	  foreach_dimension()
	    u.x[] = y > eta ? 0. : u_x(d0,d1,t+DT/2.+kWaveData.t0) * ramp_frac();
    }
    //boundary ((scalar *){u});
}

void set_wall(){
	foreach(){
		wall_mask[] = 0.0;
	}
	//boundary({wall_mask});
}


double ramp_frac () {
    return kWaveData.ramptime > 0 ? clamp(t/kWaveData.ramptime,0.,1.) : 1.;
}

double set_alpha(double x, double y, double t, double del) {

	double eta = wave(x,y,t) * ramp_frac() + kWaveData.waterlevel;
	//cout << wwelev << endl;
	if (eta < y - (del / 2.)) {
		return 0.0;
	}
	else if (eta > y + (del / 2.)) {
		return 1.0;
	}
	else {
		// Calculate volume fraction for the given cell with size del and position y
		return (eta - (y - (del / 2.))) / del;
	}
}
