#ifndef STRESS_OMEGA_TWOPHASE_H
#define STRESS_OMEGA_TWOPHASE_H
#include "rsm_common.h"
#include "diffusion.h"
#include "sander/output_htg.h"
#include "lopez/fracface.h"


scalar rhok_a[], rhok_w[], rhoe_w[], rhoe_a[];
scalar rhok[], rhoe[];

void rsm_model_debug_dump(int i, scalar* plist)
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
	foreach()
	{
		rhok[] = (alphaRho(1.0)*f[]* rhok_w[] + alphaRho(0.0)*(1.0-f[])*rhok_a[])/ (alphaRho(1.0)*f[] + alphaRho(0.0)*(1.0-f[]));
		rhoe[] = (alphaRho(1.0)*f[]* rhoe_w[] + alphaRho(0.0)*(1.0-f[])*rhoe_a[])/ (alphaRho(1.0)*f[] + alphaRho(0.0)*(1.0-f[]));
	}
	*plist = list_concat(*plist, (scalar *){rhok, rhoe, nu_t});
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
    read_double(&(kTurbConstants.k_0), "k0", tb_config);
    read_double(&(kTurbConstants.omega_0), "omega0", tb_config);

	double *dom_size = NULL;
	int rfm;
	read_num_array_1d(&(dom_size), "size", tb_config);
	read_int(&(rfm), "adapt_refine_max", tb_config);
	ls_constant.icl = rfm;
  	ls_constant.icw = (dom_size[1]-dom_size[0])/(1<<rfm);
	
	if(dom_size != NULL) {free(dom_size);}
}
#endif

event defaults(i=0) {
#if TREE
#if EMBED
  for (scalar s in {rhok, rhoe}) {
    s.restriction = restriction_embed_linear;
    s.refine = s.prolongation = refine_embed_linear;
	s.gradient = minmod2;
    s.depends = list_add (s.depends, cs);
  }
    nu_t.refine = refine_embed_linear;
    nu_t.depends = list_add (nu_t.depends, cs);
#else
	for (scalar s in {rhoe_w, rhok_w}) {
		s.refine  = refine_linear;
		s.restriction = restriction_volume_average;
		s.gradient = minmod2;
#ifdef TWO_PHASE
		s.inverse = false;
		s.depends = list_add (s.depends, f);
		f.tracers = list_add (f.tracers, s);
#endif
	}
	for (scalar s in {rhoe_a, rhok_a}) {
		s.refine  = refine_linear;
		s.restriction = restriction_volume_average;
		s.gradient = minmod2;
#ifdef TWO_PHASE
		s.inverse = true;
		s.depends = list_add (s.depends, f);
		f.tracers = list_add (f.tracers, s);
#endif
	}
	nu_t.refine  = refine_linear;

#endif
	/* On tree grids we do not directly care about the diffusivities on
refined and/or coarsend cells and faces. These should be properly
reevaluated before they appear in any computation */
	nu_t.coarsen=no_restriction;
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
        rhok_w[] = kTurbConstants.k_0 * f[];
		rhoe_w[] = kTurbConstants.omega_0 * f[];
		rhok_a[] = kTurbConstants.k_0 * (1.0 - f[]);
		rhoe_a[] = kTurbConstants.omega_0 * (1.0 - f[]);
	}
}

void bound_k(scalar rhok)
{
	foreach(){
		rhok[] = max(kTurbConstants.kMin_, rhok[]);
	}
	boundary({rhok});
}

void correct_nut() {
	bound_k(rhok_w);
	bound_k(rhok_a);
	bound(rhoe_w, kTurbConstants.omegaMin_);
	bound(rhoe_a, kTurbConstants.omegaMin_);

	foreach() 
	{
		// Blending to a wall equation
		for (scalar s in {rhok_w})
			s[] = f[] > 0.001 ?  s[] : 0.0;

		for (scalar s in {rhoe_w})
			s[] = f[] > 0.001 ?  s[] : 1e-15;

		// Blending to a wall equation
		for (scalar s in {rhok_a})
			s[] = f[] < 0.999 ?  s[] : 0.0;

		for (scalar s in {rhoe_a})
			s[] = f[] < 0.999 ?  s[] : 1e-15;
	}

	// Realizability condition
	foreach()
	{
		//double nut_w = rhok[] / (rhoe[] + 1e-15) * (1.0 - f[] * (1.0 - f[]) / 0.25) * (f[] > 0.5 ? 1.0 : 0.0);
		nu_t[] = min( (alphaRho(1.0)*f[]* rhok_w[] / (rhoe_w[] + 1e-15) + alphaRho(0.0)*(1.0-f[])*rhok_a[] / (rhoe_a[] + 1e-15))/ (alphaRho(1.0)*f[] + alphaRho(0.0)*(1.0-f[]))
		, kTurbConstants.nutMax_);
    }

	//bound2(nu_t, kTurbConstants.nutMin_, kTurbConstants.nutMax_);
	boundary({nu_t, rhok_w, rhok_a, rhoe_w, rhoe_a});

#ifdef EMBED
	foreach()
		nu_t[] *= max(cs[], 1e-14);
#endif
}

void correct_nut_bounded()
{
    bound_k(rhok_w);
	bound_k(rhok_a);
	bound(rhoe_w, kTurbConstants.omegaMin_);
	bound(rhoe_a, kTurbConstants.omegaMin_);

	foreach() 
	{
		// Blending to a wall equation
		for (scalar s in {rhok_w})
			s[] = f[] > 0.001 ?  s[] : 0.0;

		for (scalar s in {rhoe_w})
			s[] = f[] > 0.001 ?  s[] : 1e-15;

		// Blending to a wall equation
		for (scalar s in {rhok_a})
			s[] = f[] < 0.999 ?  s[] : 0.0;

		for (scalar s in {rhoe_a})
			s[] = f[] < 0.999 ?  s[] : 1e-15;
	}

	// restrict unbounded growth of tke in near potential flow.(Larsen 2018)
	// Bound oemga by stability criterion
	double a1 = 0.52;
	double lambda2 = 0.05;
	double lambda1 = 0.875;
	double Clim2 = lambda2*kTurbConstants.beta/kTurbConstants.beta_star/a1;

	foreach()
	{
	    double dudx = (u.x[1, 0, 0]-u.x[-1, 0, 0])/2.0/Delta;
	    double dudy = (u.x[0, 1, 0]-u.x[0, -1, 0])/2.0/Delta;
	    double dvdx = (u.y[1, 0, 0]-u.y[-1, 0, 0])/2.0/Delta;
	    double dvdy = (u.y[0, 1, 0]-u.y[0, -1, 0])/2.0/Delta;

		// Compute stress tensor S_{ij}
		double S_tr = (dudx + dvdy)/3.; // remove trace; divergence free
		double Sxx = dudx;
		double Syy = dvdy;
		double Sxy = 0.5*(dudy + dvdx);
		double Wxy = 0.5*(dudy - dvdx);
        double Gxx = dudx - S_tr;
        double Gyy = dvdy - S_tr;
        double Gxy = 0.5*(dudy + dvdx);

		// update |S| (to bound omega)
        double Snorm = sqrt(2*(Sxx*Sxx + Syy*Syy + 2*(Sxy*Sxy)));
		double Wnorm = sqrt(2*(2*(Wxy*Wxy)));
		double Gnorm = sqrt(2*(Gxx*Gxx + Gyy*Gyy + 2*(Gxy*Gxy)));

		// Compute p0/pw, S_{ij} and R_{ij} at current timestep are supposed to be updated.
		double p0 = 0.5 * Snorm * Snorm;
		double pw = 0.5 * Wnorm * Wnorm;
		double ts = p0/(pw+1e-15);

		double omega_star = kTurbConstants.Clim_star * Gnorm; // Wilcox 2006

		nu_t[] = min( ( alphaRho(1.0) * f[]       * rhok_w[] / (max(rhoe_w[], max(omega_star, Clim2*ts*rhoe_w[])) + 1e-15) + 
                        alphaRho(0.0) * (1.0-f[]) * rhok_a[] / (max(rhoe_a[], max(omega_star, Clim2*ts*rhoe_a[])) + 1e-15) ) / 
                        (alphaRho(1.0) * f[] + alphaRho(0.0) * (1.0-f[])), 
                        kTurbConstants.nutMax_ );

	}

	//bound2(nu_t, kTurbConstants.nutMin_, kTurbConstants.nutMax_);
	boundary({nu_t, rhok_w, rhok_a, rhoe_w, rhoe_a});

#ifdef EMBED
	foreach()
		nu_t[] *= max(cs[], 1e-14);
#endif
}


void compute_komega_mueff(face vector Dk, face vector Do, bool inverse) {
#if TREE && !AXI && FACE_FRACTION_REFINE // seems better without
    foreach_dimension() {
      Dk.x.prolongation = fraction_refine;
      Dk.x.refine = fraction_refine;
	  Do.x.prolongation = fraction_refine;
      Do.x.refine = fraction_refine;
    }
#endif

	face vector f_f[];
	face_fraction (f, f_f);

	foreach_face() {
		double mu_f = 0.5 * (alphaMuEff(f[]) + alphaMuEff(f[-1]));
		double mut_f = 0.5 * (alphaRho(f[])*nu_t[]+alphaRho(f[-1])*nu_t[-1]);
		Dk.x[] = fm.x[] * (mu_f + kTurbConstants.sigma_k*mut_f) * (inverse ? (1.0 - f_f.x[]) : f_f.x[]);
		Do.x[] = fm.x[] * (mu_f + kTurbConstants.sigma_w*mut_f) * (inverse ? (1.0 - f_f.x[]) : f_f.x[]);
	}
	boundary((scalar *){Dk, Do});
}

void compute_komega_nueff(face vector Dk, face vector Do, bool inverse) {
#if TREE && !AXI && FACE_FRACTION_REFINE // seems better without
    foreach_dimension() {
      Dk.x.prolongation = fraction_refine;
      Dk.x.refine = fraction_refine;
	  Do.x.prolongation = fraction_refine;
      Do.x.refine = fraction_refine;
    }
#endif

	face vector f_f[];
	face_fraction (f, f_f);

	foreach_face() {
		double nu_f = 0.5 * (alphaNuEff(f[]) + alphaNuEff(f[-1])) / (alphaRho(f[])+alphaRho(f[-1]));
		double nut_f = 0.5 * (alphaRho(f[])*nu_t[]+alphaRho(f[-1])*nu_t[-1]) / (alphaRho(f[])+alphaRho(f[-1]));
		Dk.x[] = fm.x[] * (nu_f + kTurbConstants.sigma_k*nut_f) * (inverse ? (1.0 - f_f.x[]) : f_f.x[]);
		Do.x[] = fm.x[] * (nu_f + kTurbConstants.sigma_w*nut_f) * (inverse ? (1.0 - f_f.x[]) : f_f.x[]);
	}
	boundary((scalar *){Dk, Do});
}

void compute_tke_srcs(
    scalar rhok, scalar rhoe,
    scalar Sk1,scalar Sp1,scalar rhov1,bool inverse
)
{
	tensor gradU[]; // grad(u)
	tensor S[], W[], tG[];

	t_rsm_biased_gradientv(u, gradU, inverse);
	
	//t_rsm_centered_gradientv(u, gradU);

	foreach(){
		// compute dui_dxi
		double tdivU = gradU.x.x[] + gradU.y.y[];

		// Compute stress tensor S_{ij}
		double S_tr = (gradU.x.x[] + gradU.y.y[])/3.; // remove trace; divergence free
		S.x.x[] = gradU.x.x[];
		S.y.y[] = gradU.y.y[];
		S.x.y[] = 0.5*(gradU.x.y[] + gradU.y.x[]);

		W.x.y[] = 0.5*(gradU.x.y[] - gradU.y.x[]);

		tG.x.x[] = S.x.x[] - S_tr;
		tG.y.y[] = S.y.y[] - S_tr;
		tG.x.y[] = S.x.y[];
		
		double alphaRho = 1.0;

#ifdef EMBED
		alphaRho *= clamp(cs[], 0.0, 1.0);
#endif

		// compute G_{ij}duidxj
		double GbyGradU = 2.0 * ( \
				tG.x.x[] * gradU.x.x[] + tG.x.y[] * gradU.x.y[] + \
				tG.x.y[] * gradU.y.x[] + tG.y.y[] * gradU.y.y[] );
			
		// Production of k, P = 2.0 \nu_t S_{ij} \frac{\partial_i}{u_j}
		//double P = alphaRho * GbyGradU / rhoeBounded2 * rhok[];
		double P = alphaRho * GbyGradU * nu_t[];

		// Correction due to divU
		//double Pcorr = -(2.0/3.0) * alphaRho * tdivU;

#ifdef TWO_PHASE
		// Buoyancy correction, Pb = - rho p_b nu_t
		// buoyancy correction, p_b = \alpha_b^* G_i drhoidxj (drhodx = (rho2-rho1)*dfdx ) (issue: unbounded growth in air side?)
		double dfdx = (f[1, 0, 0]-f[-1, 0, 0])/2.0/Delta;
		double dfdy = (f[0, 1, 0]-f[0, -1, 0])/2.0/Delta;
		double pb = kTurbConstants.alpha_b_star * (G.y * dfdy); // default : gravity goes to y
		double Pbcorr = -alphaRho * pb * nu_t[];
#else 
		double Pbcorr = 0;
#endif

		// Dissipation of k, E = \beta_\star rho k \omega
		double E = alphaRho * kTurbConstants.beta_star * rhoe[];

		// Update source terms
		Sk1[] = P + Pbcorr;
		Sp1[] = -E;
	}

#if TREE && !AXI
    for (scalar r in {rhov})
        r.prolongation = r.refine = fraction_refine;
#endif

	// Source
	foreach()
	{
        for (scalar s in {Sk1, Sp1})
            s[] *= (inverse ? (f[] < 0.999 ? 1 : 0) : (f[] > 0.001 ? 1 : 0));
	}

	boundary({Sk1, Sp1, rhov1});
}

void compute_omega_srcs(
	scalar rhok, scalar rhoe, 
	scalar Somega, scalar Spomega,
	scalar rhov, bool inverse
)
{
    tensor S[], W[], tG[];
	tensor gradU[]; // grad(u)
	vector dk[], de[];
	t_rsm_biased_gradientv(u, gradU, inverse);
	t_rsm_biased_gradient(rhok, dk, inverse);
	t_rsm_biased_gradient(rhoe, de, inverse);

	//t_rsm_centered_gradientv(u, gradU);
	//t_rsm_centered_gradient(rhok, dk);
	//t_rsm_centered_gradient(rhoe, de);

	foreach() {
		// compute dui_dxi
		double tdivU = gradU.x.x[] + gradU.y.y[] + 0.0;

		// Compute stress tensor S_{ij}
		double S_tr = (gradU.x.x[] + gradU.y.y[] + 0.0)/3.; // remove trace; divergence free
		S.x.x[] = gradU.x.x[];
		S.y.y[] = gradU.y.y[];
		//S.z.z[] = gradU.z.z[];
		S.x.y[] = 0.5*(gradU.x.y[] + gradU.y.x[]);
		S.x.z[] = 0.5*(gradU.x.z[] + gradU.z.x[]);
		//S.y.z[] = 0.5*(gradU.y.z[] + gradU.z.y[]);

		W.x.y[] = 0.5*(gradU.x.y[] - gradU.y.x[]);
		//W.x.z[] = 0.5*(gradU.x.z[] - gradU.z.x[]);
		//W.y.z[] = 0.5*(gradU.y.z[] - gradU.z.y[]);

		tG.x.x[] = S.x.x[] - S_tr;
		tG.y.y[] = S.y.y[] - S_tr;
		//tG.z.z[] = S.z.z[] - S_tr;
		tG.x.y[] = S.x.y[];
		//tG.x.z[] = S.x.z[];
		//tG.y.z[] = S.y.z[];

		// update |S| (to bound omega)
		//double Gnorm = sqrt(2*(tG.x.x[]*tG.x.x[] + tG.y.y[]*tG.y.y[] + \
		//		tG.z.z[]*tG.z.z[] + 2*(tG.x.y[]*tG.x.y[] + tG.x.z[]*tG.x.z[] + tG.y.z[]*tG.y.z[])));
		double Gnorm = sqrt(2*(tG.x.x[]*tG.x.x[] + tG.y.y[]*tG.y.y[] + 2*(tG.x.y[]*tG.x.y[])));

		double alphaRho = 1.0;

#ifdef EMBED
		alphaRho *= clamp(cs[], 0.0, 1.0);
#endif

		// compute G_{ij}duidxj
		double GbyGradU = 2.0 * ( \
				tG.x.x[] * gradU.x.x[] + tG.x.y[] * gradU.x.y[] + \
				tG.x.y[] * gradU.y.x[] + tG.y.y[] * gradU.y.y[]);

		// restrict unbounded growth of tke in near potential flow.(Larsen 2018)
		// Bound oemga by stability criterion
		double omega_star = 1e-15 + kTurbConstants.Clim_star * Gnorm; // Wilcox 2006
		double rhoeBounded = rhoe[] > omega_star ? rhoe[] : omega_star;

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

		// Update source terms
		Somega[]  = PO + PPO;
		Spomega[] =-EO;
	}

	// Source
	foreach()
	{
        for (scalar s in {Somega, Spomega})
            s[] *= (inverse ? (f[] < 0.999 ? 1 : 0) : (f[] > 0.001 ? 1 : 0));
	}

	boundary({Somega, Spomega, rhov});
}

event turbulence_correction (i++) {
	if (kVerbose>1)
		fprintf(stdout,"[%d] K-Omega: turbulence model correction.\n",i);

	// incorporate turbulence dissipation in k and omega equations
	bool mgFailed = 0;

	face vector D_kw[], D_ka[], D_ow[], D_oa[];
	scalar Sk00[], Sp00[], rhov00[];
	scalar Sk10[], Sp10[], rhov10[];
	scalar Sk01[], Sp01[], rhov01[];
	scalar Sk11[], Sp11[], rhov11[];

	compute_komega_nueff(D_kw, D_ow, false);
	compute_tke_srcs(rhok_w, rhoe_w, Sk01, Sp01, rhov01, false);
	compute_omega_srcs(rhok_w, rhoe_w, Sk00, Sp00, rhov00, false);
	boundary({rhok_w, rhoe_w});

	mgstats mgOmega_w = diffusion (rhoe_w, dt, D_ow, beta=Sp00, r=Sk00);
    mgstats mgK_w  = diffusion (rhok_w, dt, D_kw, beta=Sp01 , r=Sk01);

	compute_komega_nueff(D_ka, D_oa, true);
	compute_tke_srcs(rhok_a, rhoe_a, Sk11, Sp11, rhov11, true);
	compute_omega_srcs(rhok_a, rhoe_a, Sk10, Sp10, rhov10, true);
	boundary({rhok_a, rhoe_a});

	mgstats mgOmega_a = diffusion (rhoe_a, dt, D_oa, beta=Sp10, r=Sk10);
    mgstats mgK_a  = diffusion (rhok_a, dt, D_ka, beta=Sp11 , r=Sk11);


    if (kVerbose>0)
    {
        fprintf (stderr,
				"MG convergence for rhoe(water) reached after %d iterations, res: %g nrelax: %d\n",
				mgOmega_w.i, mgOmega_w.resa, mgOmega_w.nrelax);
        fprintf (stderr,
                "MG convergence for rhok(water) reached after %d iterations, res: %g nrelax: %d\n",
                mgK_w.i, mgK_w.resa, mgK_w.nrelax);
		fprintf (stderr,
				"MG convergence for rhoe(air)   reached after %d iterations, res: %g nrelax: %d\n",
				mgOmega_a.i, mgOmega_a.resa, mgOmega_a.nrelax);
        fprintf (stderr,
                "MG convergence for rhok(air)   reached after %d iterations, res: %g nrelax: %d\n",
                mgK_a.i, mgK_a.resa, mgK_a.nrelax);
    }

    mgFailed = mgFailed || mgOmega_w.resa > 1.0;
    mgFailed = mgFailed || mgK_w.resa > 1.0;
    mgFailed = mgFailed || mgOmega_a.resa > 1.0;
    mgFailed = mgFailed || mgK_a.resa > 1.0;

    if(mgFailed)
    {
		scalar *tmpList;
		rsm_model_output(&tmpList);
        rsm_model_debug_dump(i, tmpList);
    }

	// Compute nu_t
	// for two-phase flow merge results to rhok
	correct_nut_bounded();

	// Output stats
	if (kVerbose>0)
	{
		stats turb_nu = statsf (nu_t);
		fprintf (stderr, "[%d] RANS model nu_t  min = %e, max = %e.\n", pid(), turb_nu.min, turb_nu.max);
	}
}

// Events //
event init (i=0,last) {
	// k, omega should be initialized by user
	fprintf(stdout,"###Turbulence Model Initialization Subroutine ###\n");
	
if (!restore ("restart")) {
	set_rans_init();
	correct_nut_bounded();
	}
}

#if TREE
// mesh adaptation happens before solving viscous term
event adapt (i++) {
	if (kVerbose>1)
		fprintf(stdout,"[%d] K-omega: Correct nu_t after regridding.\n",i);
	correct_nut_bounded();
}
#endif

#endif