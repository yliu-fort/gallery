#ifndef STRESS_OMEGA_TWOPHASE_H
#define STRESS_OMEGA_TWOPHASE_H
#include "rsm_common.h"
#include "diffusion.h"
#include "sander/output_htg.h"
#include "lopez/fracface.h"

scalar Rxx_w[], Rxy_w[], Rxz_w[], Ryy_w[], Ryz_w[], Rzz_w[];
scalar Rxx_a[], Rxy_a[], Rxz_a[], Ryy_a[], Ryz_a[], Rzz_a[];
scalar rhok_a[], rhok_w[], rhoe_w[], rhoe_a[];

scalar rhok[], rhoe[];
scalar Rxx[], Rxy[], Rxz[], Ryy[], Ryz[], Rzz[];

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
		Rxx[] = (alphaRho(1.0)*f[]* Rxx_w[] + alphaRho(0.0)*(1.0-f[])*Rxx_a[])/ (alphaRho(1.0)*f[] + alphaRho(0.0)*(1.0-f[]));
		Rxy[] = (alphaRho(1.0)*f[]* Rxy_w[] + alphaRho(0.0)*(1.0-f[])*Rxy_a[])/ (alphaRho(1.0)*f[] + alphaRho(0.0)*(1.0-f[]));
		Rxz[] = (alphaRho(1.0)*f[]* Rxz_w[] + alphaRho(0.0)*(1.0-f[])*Rxz_a[])/ (alphaRho(1.0)*f[] + alphaRho(0.0)*(1.0-f[]));
		Ryy[] = (alphaRho(1.0)*f[]* Ryy_w[] + alphaRho(0.0)*(1.0-f[])*Ryy_a[])/ (alphaRho(1.0)*f[] + alphaRho(0.0)*(1.0-f[]));
		Ryz[] = (alphaRho(1.0)*f[]* Ryz_w[] + alphaRho(0.0)*(1.0-f[])*Ryz_a[])/ (alphaRho(1.0)*f[] + alphaRho(0.0)*(1.0-f[]));
		Rzz[] = (alphaRho(1.0)*f[]* Rzz_w[] + alphaRho(0.0)*(1.0-f[])*Rzz_a[])/ (alphaRho(1.0)*f[] + alphaRho(0.0)*(1.0-f[]));

		rhok[] = (alphaRho(1.0)*f[]* rhok_w[] + alphaRho(0.0)*(1.0-f[])*rhok_a[])/ (alphaRho(1.0)*f[] + alphaRho(0.0)*(1.0-f[]));
		rhoe[] = (alphaRho(1.0)*f[]* rhoe_w[] + alphaRho(0.0)*(1.0-f[])*rhoe_a[])/ (alphaRho(1.0)*f[] + alphaRho(0.0)*(1.0-f[]));
	}
	*plist = list_concat(*plist, (scalar *){rhok, rhoe, nu_t, Rxx, Rxy, Rxz, Ryy, Ryz, Rzz});
}

struct StressOmegaConstants {
    double C1;
    double alpha_v;
    double beta_v;
    double gamma_v;

    double alpha;
    double beta;
    double beta0;
	double alpha_b_star;

	double sigma_k;
	double sigma_w;
	double beta_star;
	double gamma;


	double kMin_;
	double omegaMin_;
	double nutMin_;
	double nutMax_;

	double k_0;
	double omega_0;
	double pe;
};

// C2 = (10.0/19.0)
struct StressOmegaConstants kTurbConstants = {
.C1 = 9.0/5.0,
.alpha_v = (8.0 + 10.0/19.0)/11.0,
.beta_v = (8.0*10.0/19.0 - 2.0)/11.0,
.gamma_v = (60.0*10.0/19.0 - 4.0)/55.0,

.alpha = 13.0/25.0,
.beta = 0.0708,
.beta0 = 0.0708,
.alpha_b_star = 1.36,

.sigma_k = 0.6,
.sigma_w = 0.5,
.beta_star = 0.09,

.kMin_ = 1e-15,
.omegaMin_ = 1e-10,
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

event defaults(i=0){
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
	for (scalar s in {rhoe_w, Rxx_w, Rxy_w, Rxz_w, Ryy_w, Ryz_w, Rzz_w}) {
		s.refine  = refine_linear;
		s.restriction = restriction_volume_average;
		s.gradient = minmod2;
#ifdef TWO_PHASE
		s.inverse = false;
		s.depends = list_add (s.depends, f);
		f.tracers = list_add (f.tracers, s);
#endif
	}
	for (scalar s in {rhoe_a, Rxx_a, Rxy_a, Rxz_a, Ryy_a, Ryz_a, Rzz_a}) {
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

	fprintf(stdout,"###Stress-Omega Turbulence Model ###\n");
    fprintf(stdout,"C1 = %f\n",kTurbConstants.C1 );
    fprintf(stdout,"alpha_v = %f\n",kTurbConstants.alpha_v );
    fprintf(stdout,"beta_v = %f\n",kTurbConstants.beta_v );
    fprintf(stdout,"gamma_v = %f\n",kTurbConstants.gamma_v );
    fprintf(stdout,"alpha = %f\n",kTurbConstants.alpha );
    fprintf(stdout,"beta = %f\n",kTurbConstants.beta );
    fprintf(stdout,"beta0 = %f\n",kTurbConstants.beta0 );
    fprintf(stdout,"sigma* = %f\n",kTurbConstants.sigma_k );
    fprintf(stdout,"sigma = %f\n",kTurbConstants.sigma_w );
    fprintf(stdout,"beta* = %f\n",kTurbConstants.beta_star );

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
        Rxx_w[] = -2.0/3.0*kTurbConstants.k_0*f[];
        Ryy_w[] = -2.0/3.0*kTurbConstants.k_0*f[];
        Rxy_w[] = 0.0;

		Rxz_w[] = 0.0;
		Ryz_w[] = 0.0;
		Rzz_w[] = -2.0/3.0*kTurbConstants.k_0*f[];

		rhoe_w[] = kTurbConstants.omega_0*f[];

		Rxx_a[] = -2.0/3.0*kTurbConstants.k_0 * (1 - f[]);
        Ryy_a[] = -2.0/3.0*kTurbConstants.k_0 * (1 - f[]);
        Rxy_a[] = 0.0;

		Rxz_a[] = 0.0;
		Ryz_a[] = 0.0;
		Rzz_a[] = -2.0/3.0*kTurbConstants.k_0 * (1 - f[]);

		rhoe_a[] = kTurbConstants.omega_0 * (1 - f[]);
	}
}

void bound_stresses(scalar Rxx, scalar Rxy, scalar Rxz, scalar Ryy, scalar Ryz, scalar Rzz)
{
	foreach(){
		Rxx[] = min(0, Rxx[]);
		Ryy[] = min(0, Ryy[]);
		Rzz[] = min(0, Rzz[]);

		Rxy[] = min(max(Rxy[], -sqrt(Rxx[]*Ryy[])), sqrt(Rxx[]*Ryy[]));
		#if dimension > 2
		Ryz[] = min(max(Ryz[], -sqrt(Rzz[]*Ryy[])), sqrt(Rzz[]*Ryy[]));
		Rxz[] = min(max(Rxz[], -sqrt(Rxx[]*Rzz[])), sqrt(Rxx[]*Rzz[]));
		#endif
	}
	boundary({Rxx, Rxy, Rxz, Ryy, Ryz, Rzz});
}

void correct_nut() {
	bound_stresses(Rxx_w, Rxy_w, Rxz_w, Ryy_w, Ryz_w, Rzz_w);
	bound_stresses(Rxx_a, Rxy_a, Rxz_a, Ryy_a, Ryz_a, Rzz_a);
	bound(rhoe_w, kTurbConstants.omegaMin_);
	bound(rhoe_a, kTurbConstants.omegaMin_);

	foreach() 
	{
		// Blending to a wall equation
		for (scalar s in {Rxx_w, Rxy_w, Ryy_w, Rzz_w, rhok_w})
			s[] = f[] > 0.001 ?  s[] : 0.0;

		for (scalar s in {rhoe_w})
			s[] = f[] > 0.001 ?  s[] : 1e-15;

		// Blending to a wall equation
		for (scalar s in {Rxx_a, Rxy_a, Ryy_a, Rzz_a, rhok_a})
			s[] = f[] < 0.999 ?  s[] : 0.0;

		for (scalar s in {rhoe_a})
			s[] = f[] < 0.999 ?  s[] : 1e-15;
	}

	// Realizability condition
	foreach()
	{
		rhok_w[] = max(-0.5*(Rxx_w[] + Ryy_w[] + Rzz_w[]), 0.0);
		rhok_a[] = max(-0.5*(Rxx_a[] + Ryy_a[] + Rzz_a[]), 0.0);

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

void compute_stress_srcs(
	scalar Rxx, scalar Rxy, scalar Ryy, scalar Rzz,
	scalar rhok, scalar rhoe,
	scalar Sk1, scalar Sk2, scalar Sk3, scalar Sk4,
	scalar Sp1, scalar Sp2, scalar Sp3, scalar Sp4,
	scalar rhov1, scalar rhov2, scalar rhov3, scalar rhov4,
	bool inverse
 )
{
	tensor gradU[]; // grad(u)
	tensor S[], W[], tG[];

	t_rsm_biased_gradientv(u, gradU, inverse);
	
	//t_rsm_centered_gradientv(u, gradU);

	foreach(){
		double divU = gradU.x.x[] + gradU.y.y[] + 0.0;

		// Compute stress tensor S_{ij}
		double S_tr = (gradU.x.x[] + gradU.y.y[] + 0.0)/3.; // remove trace; divergence free
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

        // Compute R production and transportation
		double Pxx, Pxy, Pyy, Pzz;
		double Yxx, Yxy, Yyy, Yzz;
		double Dxx, Dxy, Dyy, Dzz;
        Pxx = (Rxx[]*gradU.x.x[] + Rxy[]*gradU.x.y[])*2;
        Pxy =  Rxx[]*gradU.y.x[] + Rxy[]*gradU.y.y[]
			 + Rxy[]*gradU.x.x[] + Ryy[]*gradU.x.y[];

        Pyy = (Rxy[]*gradU.y.x[] + Ryy[]*gradU.y.y[])*2;

        Pzz = 0.0;

        Dxx = (Rxx[]*gradU.x.x[] + Rxy[]*gradU.y.x[])*2;
        Dxy =  Rxx[]*gradU.x.y[] + Rxy[]*gradU.y.y[]
			 + Rxy[]*gradU.x.x[] + Ryy[]*gradU.y.x[];

        Dyy = (Rxy[]*gradU.x.y[] + Ryy[]*gradU.y.y[])*2;

        Dzz = 0.0;

        double P = 0.5*(Pxx + Pyy + Pzz);

        Yxx = kTurbConstants.beta_star*kTurbConstants.C1*rhoe[]*(2.0/3.0*rhok[]) - kTurbConstants.alpha_v*(Pxx - 2.0/3.0*P)
                    - kTurbConstants.beta_v*(Dxx - 2.0/3.0*P) - kTurbConstants.gamma_v*rhok[]*(tG.x.x[]);
        Yxy =-kTurbConstants.alpha_v*(Pxy) - kTurbConstants.beta_v*(Dxy) - kTurbConstants.gamma_v*rhok[]*(S.x.y[]);

        Yyy = kTurbConstants.beta_star*kTurbConstants.C1*rhoe[]*(2.0/3.0*rhok[]) - kTurbConstants.alpha_v*(Pyy - 2.0/3.0*P)
                    - kTurbConstants.beta_v*(Dyy - 2.0/3.0*P) - kTurbConstants.gamma_v*rhok[]*(tG.y.y[]);

        Yzz = kTurbConstants.beta_star*kTurbConstants.C1*rhoe[]*(0 + 2.0/3.0*rhok[]) - kTurbConstants.alpha_v*(Pzz - 2.0/3.0*P)
                    - kTurbConstants.beta_v*(Dzz - 2.0/3.0*P);

		// Compute buoyancy correction term
		double Pbcorr_xy = 0;
		double Pbcorr_yy = 0;
#ifdef TWO_PHASE
		// buoyancy correction, p_b = \alpha_b^* G_i drhoidxj (drhodx = (rho2-rho1)*dfdx ) (issue: unbounded growth in air side?)
		double dfdx = (f[1, 0, 0]-f[-1, 0, 0])/2.0/Delta;
		double dfdy = (f[0, 1, 0]-f[0, -1, 0])/2.0/Delta;

		//Pbcorr_xy = alphaRho * kTurbConstants.alpha_b_star * G.y * dfdx * (rhok[] / (rhoe[] + 1e-15));
		//Pbcorr_yy = alphaRho * kTurbConstants.alpha_b_star * G.y * dfdy * (rhok[] / (rhoe[] + 1e-15));
#endif

		double E = 2.0/3.0*kTurbConstants.beta_star*rhoe[]*rhok[];
		double Ecorr = kTurbConstants.beta_star*kTurbConstants.C1*rhoe[];

		// Update source terms
		Sk1[] = -Pxx-Yxx + E;
        Sk2[] = -Pxy-Yxy + Pbcorr_xy;
        Sk3[] = -Pyy-Yyy + Pbcorr_yy + E;
        Sk4[] = -Pzz-Yzz + E;

		Sp1[] = -Ecorr;
        Sp2[] = -Ecorr;
        Sp3[] = -Ecorr;
        Sp4[] = -Ecorr;
	}

	#if TREE && !AXI
for (scalar r in {rhov})
    r.prolongation = r.refine = fraction_refine;
#endif

	// Source
	foreach()
	{
        for (scalar s in {Sk1, Sk2, Sk3, Sk4, Sp1, Sp2, Sp3, Sp4})
            s[] *= (inverse ? (f[] < 0.999 ? 1 : 0) : (f[] > 0.001 ? 1 : 0));
	}

	boundary({Sk1, Sk2, Sk3, Sk4, Sp1, Sp2, Sp3, Sp4, rhov1, rhov2, rhov3, rhov4});
}

void compute_omega_srcs(
	scalar Rxx, scalar Rxy, scalar Ryy, scalar Rzz,
	scalar rhok, scalar rhoe, 
	scalar Somega, scalar Spomega,
	scalar rhov, bool inverse
)
{
	tensor gradU[]; // grad(u)
	vector dk[], de[];
	t_rsm_biased_gradientv(u, gradU, inverse);
	t_rsm_biased_gradient(rhok, dk, inverse);
	t_rsm_biased_gradient(rhoe, de, inverse);

	//t_rsm_centered_gradientv(u, gradU);
	//t_rsm_centered_gradient(rhok, dk);
	//t_rsm_centered_gradient(rhoe, de);

	foreach(){

		// Compute \Chi_\omega, apply "Pope correction"
		double f_beta = 1.0;
		
		double alphaRho = 1.0;

#ifdef EMBED
		alphaRho *= clamp(cs[], 0.0, 1.0);
#endif

        // Compute the dissipation term
        // compute R_{ij}duidxj
		double RbyGradU = Rxx[] * gradU.x.x[] + Rxy[] * gradU.x.y[] + \
				          Rxy[] * gradU.y.x[] + Ryy[] * gradU.y.y[];

        double PO = kTurbConstants.alpha*RbyGradU / (rhok[] + 1e-10) * rhoe[];
        double EO = kTurbConstants.beta*f_beta*rhoe[];

		// Secondary generation of omega, PP_\omega = \sigma_d / \omega \frac{\partial_i}{k} \frac{\partial_i}{\omega}
        double dkdw = dk.x[]*de.x[]+dk.y[]*de.y[];

		double sigma_d = dkdw > 0.0 ? 0.125 : 0.0;
		double PPO = rhoe[] >= kTurbConstants.omegaMin_ ? alphaRho * sigma_d * dkdw / rhoe[] : 0.0;

		// Update source terms
		Somega[]  = PO + PPO;
		Spomega[] =  - EO;
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
		fprintf(stdout,"[%d] Stress-omega: turbulence model correction.\n",i);

	// incorporate turbulence dissipation in k and omega equations
	bool mgFailed = 0;
	scalar tmpR, tmpSU, tmpSP, tmpRHO;

	face vector D_kw[], D_ka[], D_ow[], D_oa[];
	scalar Sk00[], Sp00[], rhov00[];
	scalar Sk10[], Sp10[], rhov10[];

	scalar Sk01[], Sk02[], Sk03[], Sk04[];
	scalar Sp01[], Sp02[], Sp03[], Sp04[];
	scalar rhov01[],rhov02[],rhov03[],rhov04[];

	scalar Sk11[], Sk12[], Sk13[], Sk14[];
	scalar Sp11[], Sp12[], Sp13[], Sp14[];
	scalar rhov11[],rhov12[],rhov13[],rhov14[];

	compute_komega_nueff(D_kw, D_ow, false);
	compute_stress_srcs(
		Rxx_w, Rxy_w, Ryy_w, Rzz_w, rhok_w, rhoe_w,
		Sk01, Sk02, Sk03, Sk04,
		Sp01, Sp02, Sp03, Sp04,
		rhov01,rhov02,rhov03,rhov04, false);
	compute_omega_srcs(
		Rxx_w, Rxy_w, Ryy_w, Rzz_w, 
		rhok_w, rhoe_w, 
		Sk00, Sp00, rhov00, false);

	boundary({Rxx_w, Rxy_w, Ryy_w, Rzz_w, rhok_w, rhoe_w});

	// Diffusion
	if (kVerbose>1)
		fprintf(stdout,"[%d] Stress-Omega: diffuse turbulence variables.\n",i);

	mgstats mgOmega_w = diffusion (rhoe_w, dt, D_ow, beta=Sp00, r=Sk00);
	mgFailed = mgFailed || mgOmega_w.resa > 1.0;

	if (kVerbose>0)
	{
			fprintf (stderr,
				"MG convergence for rhoe(water) reached after %d iterations, res: %g nrelax: %d\n",
				mgOmega_w.i, mgOmega_w.resa, mgOmega_w.nrelax);
	}

    
	scalar *Rlist = (scalar *){Rxx_w, Rxy_w, Ryy_w, Rzz_w};
	scalar *SUlist = (scalar *){Sk01, Sk02, Sk03, Sk04};
	scalar *SPlist = (scalar *){Sp01, Sp02, Sp03, Sp04};
	scalar *RHOlist = (scalar *){rhov01,rhov02,rhov03,rhov04};


    for (tmpR, tmpSU, tmpSP, tmpRHO in Rlist, SUlist, SPlist, RHOlist)
    {
        mgstats mgR  = diffusion (tmpR, dt, D_kw, beta=tmpSP , r=tmpSU);

        if (kVerbose>0)
	    {
            char vxname[80]="";
			strcpy(vxname, tmpR.name);
		    const char* vname = strtok(vxname, ".");

            fprintf (stderr,
                    "MG convergence for %s reached after %d iterations, res: %g nrelax: %d\n",vname,
                    mgR.i, mgR.resa, mgR.nrelax);
        }
        mgFailed = mgFailed || mgR.resa > 1.0;
    }

	compute_komega_nueff(D_ka, D_oa, true);
	compute_stress_srcs(
		Rxx_a, Rxy_a, Ryy_a, Rzz_a, rhok_a, rhoe_a,
		Sk11, Sk12, Sk13, Sk14,
		Sp11, Sp12, Sp13, Sp14,
		rhov11,rhov12,rhov13,rhov14, true);
	compute_omega_srcs(
		Rxx_a, Rxy_a, Ryy_a, Rzz_a, 
		rhok_a, rhoe_a, 
		Sk10, Sp10, rhov10, true);

	boundary({Rxx_a, Rxy_a, Ryy_a, Rzz_a, rhok_a, rhoe_a});

	mgstats mgOmega_a = diffusion (rhoe_a, dt, D_oa, beta=Sp10, r=Sk10);
	mgFailed = mgFailed || mgOmega_a.resa > 1.0;

	if (kVerbose>0)
	{
		fprintf (stderr,
				"MG convergence for rhoe(air) reached after %d iterations, res: %g nrelax: %d\n",
				mgOmega_a.i, mgOmega_a.resa, mgOmega_a.nrelax);
	}

	scalar *Rlist2 = (scalar *){Rxx_a, Rxy_a, Ryy_a, Rzz_a};
	scalar *SUlist2 = (scalar *){Sk11, Sk12, Sk13, Sk14};
	scalar *SPlist2 = (scalar *){Sp11, Sp12, Sp13, Sp14};
	scalar *RHOlist2 = (scalar *){rhov11,rhov12,rhov13,rhov14};

    for (tmpR, tmpSU, tmpSP, tmpRHO in Rlist2, SUlist2, SPlist2, RHOlist2)
    {
        mgstats mgR  = diffusion (tmpR, dt, D_ka, beta=tmpSP , r=tmpSU);

        if (kVerbose>0)
	    {
            char vxname[80]="";
			strcpy(vxname, tmpR.name);
		    const char* vname = strtok(vxname, ".");

            fprintf (stderr,
                    "MG convergence for %s reached after %d iterations, res: %g nrelax: %d\n",vname,
                    mgR.i, mgR.resa, mgR.nrelax);
        }
        mgFailed = mgFailed || mgR.resa > 1.0;
    }


    if(mgFailed)
    {
		scalar *tmpList;
		rsm_model_output(&tmpList);
		tmpList = list_concat(tmpList, SUlist);
		tmpList = list_concat(tmpList, SPlist);
		tmpList = list_concat(tmpList, SUlist2);
		tmpList = list_concat(tmpList, SPlist2);
        rsm_model_debug_dump(i, tmpList);
    }

	// Compute nu_t
	// for two-phase flow merge results to rhok
	correct_nut();

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
	correct_nut();
	}
}

#if TREE
// mesh adaptation happens before solving viscous term
event adapt (i++) {
	if (kVerbose>1)
		fprintf(stdout,"[%d] Stress-omega: Correct nu_t after regridding.\n",i);
	correct_nut();
}
#endif

#endif