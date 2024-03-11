#ifndef KOMEGA_H
#define KOMEGA_H
#include "rsm_common.h"
#include "levelset_common.h"
#include "diffusion.h"
#include "sander/output_htg.h"

scalar rhok[], rhoe[];
scalar Snorm[],Rnorm[],Gnorm[];

mgstats mgK;
mgstats mgOmega;

void debug_dump(int i, scalar* plist)
{
	scalar vort[];
	vorticity (u, vort);

	char fname[50];
	sprintf(fname, "debug.%06d", i);
	output_htg(plist,(vector *){u}, "", fname, i, t);

	fprintf (stderr, "write output to %s\n", fname);
}

void rsm_model_output(scalar** plist)
{
	*plist = list_concat(*plist, (scalar *){rhok, rhoe, nu_t});
}

void rsm_model_output_compact(scalar** plist)
{
	*plist = list_concat(*plist, (scalar *){rhok, nu_t});
}

struct KOmegaConstants {
	double sigma_k;
	double sigma_w;
	double beta_star;
	double gamma;
	double beta;
	double Clim;
	double Clim_star;
	double alpha_b_star;
	double kMin_;
	double omegaMin_;
	double nutMin_;
	double nutMax_;

	double k_0;
	double omega_0;
	double pe;
};

struct KOmegaConstants kTurbConstants = {
.sigma_k = 0.6,
.sigma_w = 0.5,
.beta_star = 0.09,
.gamma = 13.0/25.0,
.beta = 0.0708,
.Clim = 7.0/8.0,
.Clim_star = 7.0/8.0/0.3,
.alpha_b_star = 1.36,
.kMin_ = 1e-15,
.omegaMin_ = 1e-10, // In order to make mg solver happy, use a larger sup-bound for omega...
.nutMin_ = 1e-15,
.nutMax_ = 1e3,

.k_0 = 1e-4,
.omega_0 = 420.0,
.pe = 1e10,
};

#ifdef JSON_UTILS_H
void refresh_rsm_model_by_json(struct rsm_JSONFile_t *jf)
{
    cJSON *tb_config = jf->json;

	// Update k_0
    read_double(&(kTurbConstants.k_0), "k0", tb_config);

	// Update \omega_0
    read_double(&(kTurbConstants.omega_0), "omega0", tb_config);

	// Update interface cell level and width
	double *dom_size = NULL;
	int rfm;
	read_num_array_1d(&(dom_size), "size", tb_config);
	read_int(&(rfm), "adapt_refine_max", tb_config);
	ls_constant.icl = rfm;
	if(dom_size){
		ls_constant.icw = (dom_size[1]-dom_size[0])/(1<<rfm);
	}else{
		fprintf(stderr, "[RSM Model]:Failed to load domain size! Please check the system json file.\n");
		ls_constant.icw = 0;
	}
  	
	if(dom_size != NULL) {free(dom_size);}
}
#endif

event defaults(i=0){
#if TREE
#if EMBED
  for (scalar s in {rhok, rhoe}) {

    s.restriction = restriction_embed_linear;
    s.refine = s.prolongation = refine_embed_linear;
    s.depends = list_add (s.depends, cs);
	s.gradient = minmod;
  }
    nu_t.refine = refine_embed_linear;
    nu_t.depends = list_add (nu_t.depends, cs);
	nu_t.coarsen = no_restriction;
#else
	rhok.refine  = refine_linear;
	rhoe.refine  = refine_linear;

	rhok.restriction = restriction_volume_average;
	rhoe.restriction = restriction_volume_average;

	rhok.gradient = minmod2;
	rhoe.gradient = minmod2;

	nu_t.refine  = refine_linear;
	nu_t.coarsen = nu_t.restriction = restriction_volume_average;
#endif
	/* On tree grids we do not directly care about the diffusivities on
refined and/or coarsend cells and faces. These should be properly
reevaluated before they appear in any computation */
	foreach_dimension(){
		mu_t.x.coarsen=no_restriction;
	}
#endif

	fprintf(stdout,"###K-Omega Turbulence Model ###\n");
	fprintf(stdout,"sigmaK = %f\n",kTurbConstants.sigma_k);
	fprintf(stdout,"sigmaW = %f\n",kTurbConstants.sigma_w);
	fprintf(stdout,"beta* = %f\n",kTurbConstants.beta_star);
	fprintf(stdout,"gamma = %f\n",kTurbConstants.gamma);
	fprintf(stdout,"beta = %f\n",kTurbConstants.beta);
	fprintf(stdout,"Clim = %f\n",kTurbConstants.Clim);
	fprintf(stdout,"alpha_beta* = %f\n",kTurbConstants.alpha_b_star);

	fprintf(stdout,"kMin = %g\n",kTurbConstants.kMin_);
	fprintf(stdout,"omegaMin = %g\n",kTurbConstants.omegaMin_);
	fprintf(stdout,"nutMin = %g\n",kTurbConstants.nutMin_);
	fprintf(stdout,"nutMax = %g\n",kTurbConstants.nutMax_);
	fprintf(stdout,"k0 = %g\n",kTurbConstants.k_0);
	fprintf(stdout,"omega0 = %g\n",kTurbConstants.omega_0);
	fprintf(stdout,"cell peclet number = %g\n",kTurbConstants.pe);
	fprintf(stdout,"###############################\n");
}

void set_rans_init()
{
	// rans model
	foreach() {
		rhok[] = kTurbConstants.k_0;
		rhoe[] = kTurbConstants.omega_0;
		nu_t[] = kTurbConstants.k_0/kTurbConstants.omega_0;
	}

	//boundary((scalar*){rhok, rhoe, nu_t});
}

event init (i=0,last) {
	// k, omega should be initialized by user
	fprintf(stdout,"###Turbulence Model Initialization Subroutine ###\n");
if (!restore ("restart")) {
	set_rans_init();
	}
}

/*
rhok[right] = neumann (0);
rhok[left]  = neumann (0);
rhoe[right] = neumann (0);
rhoe[left]  = neumann (0);

#if AXI
rhok[top]    = neumann (0);
rhoe[top]    = neumann (0);
#else // !AXI
#  if dimension > 1
rhok[top]    = neumann (0);
rhok[bottom] = neumann (0);
rhoe[top]    = neumann (0);
rhoe[bottom] = neumann (0);
#  endif
#  if dimension > 2
rhok[front]  = neumann (0);
rhok[back]   = neumann (0);
rhoe[front]  = neumann (0);
rhoe[back]   = neumann (0);
#  endif
#endif // !AXI
*/

void smooth_nut()
{
  /**
  When using smearing of the density jump, we initialise *sf* with the
  vertex-average of *f*. */
	scalar nuTilda[];
	//boundary({nu_t});

	#if dimension <= 2
	  foreach()
	    nuTilda[] = (4.*nu_t[] + 
		    2.*(nu_t[0,1] + nu_t[0,-1] + nu_t[1,0] + nu_t[-1,0]) +
		    nu_t[-1,-1] + nu_t[1,-1] + nu_t[1,1] + nu_t[-1,1])/16.;
	#else // dimension == 3
	  foreach()
	    nuTilda[] = (8.*nu_t[] +
		    4.*(nu_t[-1] + nu_t[1] + nu_t[0,1] + nu_t[0,-1] + nu_t[0,0,1] + nu_t[0,0,-1]) +
		    2.*(nu_t[-1,1] + nu_t[-1,0,1] + nu_t[-1,0,-1] + nu_t[-1,-1] + 
			nu_t[0,1,1] + nu_t[0,1,-1] + nu_t[0,-1,1] + nu_t[0,-1,-1] +
			nu_t[1,1] + nu_t[1,0,1] + nu_t[1,-1] + nu_t[1,0,-1]) +
		    nu_t[1,-1,1] + nu_t[-1,1,1] + nu_t[-1,1,-1] + nu_t[1,1,1] +
		    nu_t[1,1,-1] + nu_t[-1,-1,-1] + nu_t[1,-1,-1] + nu_t[-1,-1,1])/64.;
	#endif

	foreach()
		nu_t[] = nuTilda[];

	//boundary({nu_t});
}

void bound_komega()
{
	bound(rhok, kTurbConstants.kMin_);
	bound(rhoe, kTurbConstants.omegaMin_);

#ifdef EMBED
	foreach()
	{
	  rhok[] = cs[] > 0.01 ? rhok[] : 0.0;
	  rhoe[] = cs[] > 0.01 ? rhoe[] : kTurbConstants.omegaMin_;
	}
#endif

	boundary({rhok, rhoe});
}

void correct_nut()
{
	foreach() {
		nu_t[] = rhok[]/(rhoe[] + 1e-15);
    }
	bound2(nu_t, kTurbConstants.nutMin_, kTurbConstants.nutMax_);
	//boundary({nu_t});

#ifdef EMBED
	foreach()
		nu_t[] *= max(cs[], 1e-14);
#endif
}

#ifndef TWO_PHASE
void correct_nut_bounded()
{
	foreach()
	{
		// Wilcox 2006
		double omega_star = kTurbConstants.Clim_star * Gnorm[];
		double rhoeBounded = rhoe[] > omega_star ? rhoe[] : omega_star;

		nu_t[] = rhok[]/(rhoeBounded+1e-15);

	}
	bound2(nu_t, kTurbConstants.nutMin_, kTurbConstants.nutMax_);

#ifdef EMBED
	foreach()
		nu_t[] *= max(cs[], 1e-14);
#endif

	//boundary({nu_t});
}
#else
void correct_nut_bounded()
{
	// restrict unbounded growth of tke in near potential flow.(Larsen 2018)
	// Bound oemga by stability criterion
	double a1 = 0.52;
	double lambda2 = 0.05;
	double lambda1 = 0.875;
	double Clim2 = lambda2*kTurbConstants.beta/kTurbConstants.beta_star/a1;

	tensor gradU[];
	tensor S[],R[],tG[];
	t_rsm_centered_gradientv(u, gradU);

	foreach()
	{
		double S_tr = (gradU.x.x[] + gradU.y.y[])/3.; // remove trace; divergence free
		S.x.x[] = gradU.x.x[];
		S.y.y[] = gradU.y.y[];
		S.x.y[] = 0.5*(gradU.x.y[] + gradU.y.x[]);

		R.x.y[] = 0.5*(gradU.x.y[] - gradU.y.x[]);

		tG.x.x[] = S.x.x[] - S_tr;
		tG.y.y[] = S.y.y[] - S_tr;
		tG.x.y[] = S.x.y[];

		// update |S| (to bound omega)
		Snorm[] = sqrt(2*(S.x.x[]*S.x.x[] + S.y.y[]*S.y.y[] + 2*(S.x.y[]*S.x.y[])));
		Rnorm[] = sqrt(2*(2*(R.x.y[]*R.x.y[])));
		Gnorm[] = sqrt(2*(tG.x.x[]*tG.x.x[] + tG.y.y[]*tG.y.y[] + 2*(tG.x.y[]*tG.x.y[])));

		// Compute p0/pw, S_{ij} and R_{ij} at current timestep are supposed to be updated.
		double p0 = 0.5 * Snorm[] * Snorm[];
		double pw = 0.5 * Rnorm[] * Rnorm[];
		double ts = p0/(pw+1e-15);

		double omega_star = kTurbConstants.Clim_star * Gnorm[]; // Wilcox 2006
		double omega_star2 = Clim2*ts*rhoe[];

		double rhoeBounded = rhoe[] > omega_star ? rhoe[] : omega_star;
		rhoeBounded = rhoeBounded > omega_star2 ? rhoeBounded : omega_star2;

		nu_t[] = rhok[]/(rhoeBounded+1e-15);

	}

	bound2(nu_t, kTurbConstants.nutMin_, kTurbConstants.nutMax_);

#ifdef EMBED
	foreach()
		nu_t[] *= max(cs[], 1e-14);
#endif

	//boundary({nu_t});
}
#endif


#if dimension == 2
void compute_komega_srcs(scalar Sk, scalar Somega, scalar Spk, scalar Spomega)
{
	tensor gradU[]; // grad(u)
	tensor S[];
	tensor R[];
	tensor tG[];
	vector df[], dk[], de[];

	t_rsm_centered_gradientv(u, gradU);
	t_rsm_centered_gradient(f, df);
	t_rsm_centered_gradient(rhok, dk);
	t_rsm_centered_gradient(rhoe, de);

	// Compute gradU and S_{ij}
	foreach(){
		// Compute stress tensor S_{ij}
		double S_tr = (gradU.x.x[] + gradU.y.y[])/3.; // remove trace; divergence free
		S.x.x[] = gradU.x.x[];
		S.y.y[] = gradU.y.y[];
		S.x.y[] = 0.5*(gradU.x.y[] + gradU.y.x[]);

		R.x.y[] = 0.5*(gradU.x.y[] - gradU.y.x[]);

		tG.x.x[] = S.x.x[] - S_tr;
		tG.y.y[] = S.y.y[] - S_tr;
		tG.x.y[] = S.x.y[];

		// update |S| (to bound omega)
		Snorm[] = sqrt(2*(S.x.x[]*S.x.x[] + S.y.y[]*S.y.y[] + 2*(S.x.y[]*S.x.y[])));
		Rnorm[] = sqrt(2*(2*(R.x.y[]*R.x.y[])));
		Gnorm[] = sqrt(2*(tG.x.x[]*tG.x.x[] + tG.y.y[]*tG.y.y[] + 2*(tG.x.y[]*tG.x.y[])));
		
		double alphaRho = 1.0;

#ifdef EMBED
		alphaRho *= clamp(cs[], 0.0, 1.0);
#endif

		// compute G_{ij}duidxj
		double GbyGradU = 2.0 * ( \
				tG.x.x[] * gradU.x.x[] + tG.x.y[] * gradU.x.y[] + \
				tG.x.y[] * gradU.y.x[] + tG.y.y[] * gradU.y.y[] );
			
		// restrict unbounded growth of tke in near potential flow.(Larsen 2018)
		// Bound oemga by stability criterion
		double omega_star = 1e-15 + kTurbConstants.Clim_star * Gnorm[]; // Wilcox 2006
		double rhoeBounded = rhoe[] > omega_star ? rhoe[] : omega_star;

		// Production of k, P = 2.0 \nu_t S_{ij} \frac{\partial_i}{u_j}
		//double P = alphaRho * GbyGradU / rhoeBounded2 * rhok[];
		double P = alphaRho * GbyGradU * nu_t[];

		// Correction due to divU (from OpenFOAM, don't really understand why)
		//double Pcorr = -(2.0/3.0) * alphaRho * tdivU;

#if defined (TWO_PHASE) && defined (RSM_BUOYANCY)
		// Buoyancy correction, Pb = - rho p_b nu_t
		// buoyancy correction, p_b = \alpha_b^* G_i drhoidxj (drhodx = (rho2-rho1)*dfdx ) (issue: unbounded growth in air side?)
		double pb = kTurbConstants.alpha_b_star * (G.y*df.y[]); // default : gravity goes to y
		double Pbcorr = -alphaRho * pb * nu_t[];
#else 
		double Pbcorr = 0;
#endif

		// Dissipation of k, E = \beta_\star rho k \omega
		double E  = alphaRho * kTurbConstants.beta_star * rhoe[];

		// Production of omega, P_\omega = \gamma \omega / \k  P
		double PO = alphaRho * kTurbConstants.gamma * (rhoe[] / rhoeBounded) *  GbyGradU;

		// Dissipation of omega, E_\omega = \beta rho \omega^2
		double EO = alphaRho * kTurbConstants.beta * rhoe[];

		// Secondary generation of omega, PP_\omega = \sigma_d / \omega \frac{\partial_i}{k} \frac{\partial_i}{\omega}
		// compute dki_dxi domegai_dxi
		// notice the treatment to \nabla\omega/\omega (assume \omega is bounded to a positive value properly) to avoid a explicit division to \omega which could lead to numerical blow-up.
		double dkdw = dk.x[]*de.x[]+dk.y[]*de.y[];
		double sigma_d = dkdw > 0.0? 0.125 : 0.0;
		double PPO = rhoe[] > kTurbConstants.omegaMin_ ? alphaRho * sigma_d * dkdw / rhoe[] : 0.0;

		// Interface damping term (Ergorov interface damping scheme)
#if defined (TWO_PHASE) && defined (RSM_ERGOROV_INTERFACE_DAMPING) 
		double inv_beta_w2 = 1.0/kTurbConstants.beta_star/sq(kTurbConstants.icw);
		double wInt = 500*6*alphaNuEff(f[])*inv_beta_w2; // B = 500 (Gada et al. 2017)
		PO += dist[] * alphaRho * kTurbConstants.beta * wInt * wInt;
#endif

		// Update source terms
		Sk[] = P + Pbcorr;
		Spk[]=-E;
		Somega[]  = PO + PPO;
		Spomega[] =-EO;

	}
	//boundary({Sk, Somega, Spk, Spomega});
}
#else
void compute_komega_srcs(scalar Sk, scalar Somega, scalar Spk, scalar Spomega)
{
	tensor gradU[]; // grad(u)
	symmetric tensor S[];
	tensor R[];
	symmetric tensor tG[];
	vector df[], dk[], de[];

	// compute dui_dxj
	t_rsm_centered_gradientv(u, gradU);
	t_rsm_centered_gradient(f, df);
	t_rsm_centered_gradient(rhok, dk);
	t_rsm_centered_gradient(rhoe, de);

	// interface damping constant
	//double inv_beta_w2 = 1.0/kTurbConstants.beta_star/sq(kTurbConstants.icw);
	//compute_distance_to_interface();

	// Compute gradU and S_{ij}
	foreach(){
		// Compute stress tensor S_{ij}
		double S_tr = (gradU.x.x[] + gradU.y.y[] + gradU.z.z[])/3.; // remove trace; divergence free
		S.x.x[] = gradU.x.x[];
		S.y.y[] = gradU.y.y[];
		S.z.z[] = gradU.z.z[];
		S.x.y[] = 0.5*(gradU.x.y[] + gradU.y.x[]);
		S.x.z[] = 0.5*(gradU.x.z[] + gradU.z.x[]);
		S.y.z[] = 0.5*(gradU.y.z[] + gradU.z.y[]);

		R.x.y[] = 0.5*(gradU.x.y[] - gradU.y.x[]);
		R.x.z[] = 0.5*(gradU.x.z[] - gradU.z.x[]);
		R.y.z[] = 0.5*(gradU.y.z[] - gradU.z.y[]);

		tG.x.x[] = S.x.x[] - S_tr;
		tG.y.y[] = S.y.y[] - S_tr;
		tG.z.z[] = S.z.z[] - S_tr;
		tG.x.y[] = S.x.y[];
		tG.x.z[] = S.x.z[];
		tG.y.z[] = S.y.z[];

		// update |S| (to bound omega)
		Snorm[] = sqrt(2*(S.x.x[]*S.x.x[] + S.y.y[]*S.y.y[] + \
				S.z.z[]*S.z.z[] + 2*(S.x.y[]*S.x.y[] + S.x.z[]*S.x.z[] + S.y.z[]*S.y.z[])));

		Rnorm[] = sqrt(2*(2*(R.x.y[]*R.x.y[] + R.x.z[]*R.x.z[] + R.y.z[]*R.y.z[])));

		Gnorm[] = sqrt(2*(tG.x.x[]*tG.x.x[] + tG.y.y[]*tG.y.y[] + \
				tG.z.z[]*tG.z.z[] + 2*(tG.x.y[]*tG.x.y[] + tG.x.z[]*tG.x.z[] + tG.y.z[]*tG.y.z[])));
		
		double alphaRho = 1.0;

#ifdef EMBED
		alphaRho *= clamp(cs[], 0.0, 1.0);
#endif

		// compute G_{ij}duidxj
		double GbyGradU = 2.0 * ( \
				tG.x.x[] * gradU.x.x[] + tG.x.y[] * gradU.x.y[] + tG.x.z[] * gradU.x.z[] + \
				tG.x.y[] * gradU.y.x[] + tG.y.y[] * gradU.y.y[] + tG.y.z[] * gradU.y.z[] + \
				tG.x.z[] * gradU.z.x[] + tG.y.z[] * gradU.z.y[] + tG.z.z[] * gradU.z.z[]);
			
		// restrict unbounded growth of tke in near potential flow.(Larsen 2018)
		// Bound oemga by stability criterion
		double omega_star = 1e-15 + kTurbConstants.Clim_star * Gnorm[]; // Wilcox 2006
		double rhoeBounded = rhoe[] > omega_star ? rhoe[] : omega_star;

		// Production of k, P = 2.0 \nu_t S_{ij} \frac{\partial_i}{u_j}
		//double P = alphaRho * GbyGradU / rhoeBounded2 * rhok[];
		double P = alphaRho * GbyGradU * nu_t[];

		// Correction due to divU
		//double Pcorr = -(2.0/3.0) * alphaRho * tdivU;

#if defined (TWO_PHASE) && defined (RSM_BUOYANCY)
		// Buoyancy correction, Pb = - rho p_b nu_t
		// buoyancy correction, p_b = \alpha_b^* G_i drhoidxj (drhodx = (rho2-rho1)*dfdx ) (issue: unbounded growth in air side?)
		double pb = kTurbConstants.alpha_b_star * (G.y*df.y[]); // default : gravity goes to y
		double Pbcorr = -alphaRho * pb * nu_t[];
#else 
		double Pbcorr = 0;
#endif

		// Dissipation of k, E = \beta_\star rho k \omega
		double E  = alphaRho * kTurbConstants.beta_star * rhoe[];

		// Production of omega, P_\omega = \gamma \omega / \k  P
		double PO = alphaRho * kTurbConstants.gamma * (rhoe[] / rhoeBounded) *  GbyGradU;
		
		// Dissipation of omega, E_\omega = \beta rho \omega^2
		double EO = alphaRho * kTurbConstants.beta * rhoe[];

		// Secondary generation of omega, PP_\omega = \sigma_d / \omega \frac{\partial_i}{k} \frac{\partial_i}{\omega}
		// compute dki_dxi domegai_dxi
		// notice the treatment to \nabla\omega/\omega (assume \omega is bounded to a positive value properly) to avoid a explicit division to \omega which could lead to numerical blow-up.
		double dkdw = dk.x[]*de.x[]+dk.y[]*de.y[]+dk.z[]*de.z[];
		double sigma_d = dkdw > 0.0? 0.125 : 0.0;
		double PPO = rhoe[] > kTurbConstants.omegaMin_ ? alphaRho * sigma_d * dkdw / rhoe[] : 0.0;

		// Interface damping term (Ergorov interface damping scheme)
#if defined (TWO_PHASE) && defined (RSM_ERGOROV_INTERFACE_DAMPING) 
		double inv_beta_w2 = 1.0/kTurbConstants.beta_star/sq(kTurbConstants.icw);
		double wInt = 500*6*alphaNuEff(f[])*inv_beta_w2; // B = 500 (Gada et al. 2017)
		PO += dist[] * alphaRho * kTurbConstants.beta * wInt * wInt;
#endif

		// Update source terms
		Sk[] = P + Pbcorr;
		Spk[]=-E;
		Somega[]  = PO + PPO;
		Spomega[] =-EO;

	}
	//boundary({Sk, Somega, Spk, Spomega});
}
#endif
void compute_komega_mueff(face vector Dk, face vector Domega)
{
	foreach_face() {
		double mu_f = 0.5 * (alphaMuEff(f[]) + alphaMuEff(f[-1]));
		double mut_f = 0.5 * (alphaRho(f[])*nu_t[]+alphaRho(f[-1])*nu_t[-1]);
		Dk.x[] = fm.x[] * (mu_f + kTurbConstants.sigma_k*mut_f);
		Domega.x[] = fm.x[] * (mu_f + kTurbConstants.sigma_w*mut_f);
	}
	//boundary({Dk,Domega});
}

event turbulence_correction (i++) {
	if (kVerbose>1)
		fprintf(stdout,"[%d] K-Omega: turbulence model correction.\n",i);

	// incorporate turbulence dissipation in k and omega equations
	face vector D_k[], D_omega[];
	compute_komega_mueff(D_k, D_omega);
	scalar src_k[], src_omega[], sp_k[], sp_omega[];
	compute_komega_srcs(src_k, src_omega, sp_k, sp_omega);
	scalar rhov1[], rhov2[];

	// Advection
	if (kVerbose>1)
		fprintf(stdout,"[%d] K-Omega: advect turbulence variables.\n",i);
	advection ((scalar *){rhok, rhoe}, uf, dt);
	bound_komega();
	
	// Source
	foreach()
	{
		rhov1[] = cm[]*alphaRho(f[]);
		rhov2[] = cm[]*alphaRho(f[]);
		src_omega[] *= rhov1[];
		sp_omega[] *= rhov1[];
		src_k[] *= rhov2[];
		sp_k[] *= rhov2[];
#if EMBED
		rhov1[] *= max(cs[], 1e-14);
		rhov2[] *= max(cs[], 1e-14);
#endif
	}
	//boundary({rhok, rhoe});

	// Diffusion
	if (kVerbose>1)
		fprintf(stdout,"[%d] K-Omega: diffuse turbulence variables.\n",i);

	mgOmega = diffusion (rhoe, dt, D_omega, beta=sp_omega, r=src_omega, theta=rhov1);
	mgK     = diffusion (rhok, dt, D_k    , beta=sp_k    , r=src_k    , theta=rhov2);

	if (kVerbose>0)
	{
		fprintf (stderr,
				"MG convergence for rhok reached after %d iterations, res: %g nrelax: %d\n",
				mgK.i, mgK.resa, mgK.nrelax);
		fprintf (stderr,
				"MG convergence for rhoe reached after %d iterations, res: %g nrelax: %d\n",
				mgOmega.i, mgOmega.resa, mgOmega.nrelax);

		if(mgK.resa > 1.0 || mgOmega.resa > 1.0)
		{
			//set_rans_init();
		  debug_dump(i,(scalar *){f, nu_t,rhok, rhoe, src_k, src_omega});
		}
	}

	// Compute nu_t
	bound_komega();
	correct_nut_bounded();

	// Output stats
	if (kVerbose>0)
	{
		stats turb_k = statsf (rhok);
		stats turb_o = statsf (rhoe);
		stats turb_nu = statsf (nu_t);
		fprintf (stderr, "[%d] RANS model k     min = %e, max = %e.\n", pid(), turb_k.min, turb_k.max);
		fprintf (stderr, "[%d] RANS model omega min = %e, max = %e.\n", pid(), turb_o.min, turb_o.max);
		fprintf (stderr, "[%d] RANS model nu_t  min = %e, max = %e.\n", pid(), turb_nu.min, turb_nu.max);
	}
}

#if TREE
// mesh adaptation happens before solving viscous term
event adapt (i++) {
	if (kVerbose>1)
		fprintf(stdout,"[%d] K-Omega: Correct nu_t after regridding.\n",i);
	//boundary({rhok, rhoe});
	correct_nut_bounded();
}
#endif

// overload viscosity eventcorrect_nut
event viscous_term(i++) {
	//if (kVerbose>1)
	//	fprintf(stdout,"[%d] RSM model: apply turbulent viscosity.\n",i);

	// add nu_t to viscosity term
	// (ISSUE) for two phase fluid viscosity at the interface cell has to be set to a small constant
	// or the poisson solver becomes unstable.
	foreach_face() {
		double mu_f = 0.5 * (alphaMu(f[]) + alphaMu(f[-1]));
#ifndef RSM_NO_INTERFACE_DIFFUSION
	double mut_f = 0.5 * (alphaRho(f[])*nu_t[]+alphaRho(f[-1])*nu_t[-1]);
#else
	double mut_f = 0.5 * (alphaRho(f[])*nu_t[]*(1.0 - dist[])+alphaRho(f[-1])*nu_t[-1]*(1.0 - dist[-1]));
#endif
		mu_t.x[] = mu_f + mut_f;  // Face center.  Multiplied by rho to convert to dynamic viscosity
	}
	//boundary((scalar *){mu_t});

	// add to dynamic viscosity
	face vector muv = mu;
	foreach_face()
	  muv.x[] = fm.x[] * mu_t.x[];
}
#endif