#ifndef STRESS_OMEGA_H
#define STRESS_OMEGA_H
#include "rsm_common.h"
#include "diffusion.h"
#include "sander/output_htg.h"

scalar rhok[], rhoe[];
scalar Rxx[], Rxy[], Rxz[], Ryy[], Ryz[], Rzz[];

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
	*plist = list_concat(*plist, (scalar *){rhok, rhoe, nu_t, Rxx, Rxy, Rxz, Ryy, Ryz, Rzz});
}

void rsm_model_output_compact(scalar** plist)
{
	*plist = list_concat(*plist, (scalar *){rhok, nu_t});
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

#define C2 (10.0/19.0)
struct StressOmegaConstants kTurbConstants = {
.C1 = 9.0/5.0,
.alpha_v = (8.0 + C2)/11.0,
.beta_v = (8.0*C2 - 2.0)/11.0,
.gamma_v = (60.0*C2 - 4.0)/55.0,

.alpha = 13.0/25.0,
.beta = 0.0708, // (TODO)
.beta0 = 0.0708,
.alpha_b_star = 1.36,

.sigma_k = 0.6,
.sigma_w = 0.5,
.beta_star = 0.09,

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
  for (scalar s in {rhok, rhoe, Rxx, Rxy, Rxz, Ryy, Ryz, Rzz}) {
    s.restriction = restriction_embed_linear;
    s.refine = s.prolongation = refine_embed_linear;
	s.gradient = minmod2;
    s.depends = list_add (s.depends, cs);
  }
    nu_t.refine = refine_embed_linear;
    nu_t.depends = list_add (nu_t.depends, cs);
#else
	for (scalar s in {rhok, rhoe, Rxx, Rxy, Rxz, Ryy, Ryz, Rzz}) {
		s.refine  = refine_linear;
		s.restriction = restriction_volume_average;
		s.gradient = minmod2;
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
	fprintf(stdout,"interface cell refinement level = %g\n",ls_constant.icl);
	fprintf(stdout,"interface cell width = %g\n",ls_constant.icw);
	fprintf(stdout,"###############################\n");
}

void set_rans_init()
{
	// rans model
	foreach() {
        Rxx[] = -2.0/3.0*kTurbConstants.k_0;
        Ryy[] = -2.0/3.0*kTurbConstants.k_0;
        Rxy[] = 0.0;

		Rxz[] = 0.0;
		Ryz[] = 0.0;
		Rzz[] = -2.0/3.0*kTurbConstants.k_0;

		rhoe[] = kTurbConstants.omega_0;
	}
	//boundary((scalar*){Rxx, Rxy, Rxz, Ryy, Ryz, Rzz, rhoe});
}

event init (i=0,last) {
	// k, omega should be initialized by user
	fprintf(stdout,"###Turbulence Model Initialization Subroutine ###\n");
if (!restore ("restart")) {
	set_rans_init();
	correct_nut();
	}
}

void correct_nut()
{
	bound(rhoe, kTurbConstants.omegaMin_);
	bound2(Rxx, -1e10,0);
	bound2(Ryy, -1e10,0);
	bound2(Rzz, -1e10,0);

	// Realizability condition
	foreach()
	{
		Rxy[] = min(max(Rxy[], -sqrt(Rxx[]*Ryy[])), sqrt(Rxx[]*Ryy[]));
		#if dimension > 2
		Ryz[] = min(max(Ryz[], -sqrt(Rzz[]*Ryy[])), sqrt(Rzz[]*Ryy[]));
		Rxz[] = min(max(Rxz[], -sqrt(Rxx[]*Rzz[])), sqrt(Rxx[]*Rzz[]));
		#endif
	}

	foreach() {
#ifdef EMBED
		Rxx[] = cs[] > 0.01 ? Rxx[] : 0.0;
		Rxy[] = cs[] > 0.01 ? Rxy[] : 0.0;
		Rxz[] = cs[] > 0.01 ? Rxz[] : 0.0;
		Ryy[] = cs[] > 0.01 ? Ryy[] : 0.0;
		Ryz[] = cs[] > 0.01 ? Ryz[] : 0.0;
		Rzz[] = cs[] > 0.01 ? Rzz[] : 0.0;
		rhoe[] = cs[] > 0.01 ? rhoe[] : kTurbConstants.omegaMin_;
#endif
		rhok[] = -0.5*(Rxx[] + Ryy[] + Rzz[]);
		nu_t[] = rhok[]/(rhoe[] + 1e-15);
    }

	bound2(nu_t, kTurbConstants.nutMin_, kTurbConstants.nutMax_);
	//boundary({rhok, rhoe, nu_t});

#ifdef EMBED
	foreach()
		nu_t[] *= max(cs[], 1e-14);
#endif
}

#if dimension == 2
void compute_stress_omega_srcs(
	scalar Sk01, scalar Sk02, scalar Sk03,
	scalar Sk04,
	scalar Spk01, scalar Spk02, scalar Spk03,
	scalar Spk04,
scalar Somega, scalar Spomega
)
{
	tensor gradU[]; // grad(u)
	tensor S[], W[], tG[];
	vector df[], dk[], de[];

	t_rsm_centered_gradientv(u, gradU);
	t_rsm_centered_gradient(rhok, dk);
	t_rsm_centered_gradient(rhoe, de);

	foreach(){
		// compute dui_dxi
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

#ifdef TWO_PHASE
		// buoyancy correction, p_b = \alpha_b^* G_i drhoidxj (drhodx = (rho2-rho1)*dfdx ) (issue: unbounded growth in air side?)
		double dfdx = (f[1, 0, 0]-f[-1, 0, 0])/2.0/Delta;
		double dfdy = (f[0, 1, 0]-f[0, -1, 0])/2.0/Delta;

		double Pbcorr_xy = alphaRho * kTurbConstants.alpha_b_star * G.y * dfdx * nu_t[];
		double Pbcorr_yy = alphaRho * kTurbConstants.alpha_b_star * G.y * dfdy * nu_t[];
#else 
		double Pbcorr_xy = 0;
		double Pbcorr_yy = 0;
#endif

		double E = 2.0/3.0*kTurbConstants.beta_star*rhoe[]*rhok[];
		double Ecorr = kTurbConstants.beta_star*kTurbConstants.C1*rhoe[];

        // Compute the dissipation term
        // compute R_{ij}duidxj
		double RbyGradU = Rxx[] * gradU.x.x[] + Rxy[] * gradU.x.y[] + \
				          Rxy[] * gradU.y.x[] + Ryy[] * gradU.y.y[];

        double PO = kTurbConstants.alpha*RbyGradU / (rhok[] + 1e-15) * rhoe[];
        double EO = kTurbConstants.beta*rhoe[];

		// Secondary generation of omega, PP_\omega = \sigma_d / \omega \frac{\partial_i}{k} \frac{\partial_i}{\omega}
		// compute dki_dxi domegai_dxi
		// notice the treatment to \nabla\omega/\omega (assume \omega is bounded to a positive value properly) to avoid a explicit division to \omega which could lead to numerical blow-up.
        double dkdw = dk.x[]*de.x[]+dk.y[]*de.y[];
		double sigma_d = dkdw > 0.0? 0.125 : 0.0;
		double PPO = rhoe[] > kTurbConstants.omegaMin_ ? alphaRho * sigma_d * dkdw / rhoe[] : 0.0;

		// Update source terms
		Sk01[] = -Pxx-Yxx + E;
        Sk02[] = -Pxy-Yxy + Pbcorr_xy;
        Sk03[] = -Pyy-Yyy + E + Pbcorr_yy;
        Sk04[] = -Pzz-Yzz + E;

		Spk01[] = -Ecorr;
        Spk02[] = -Ecorr;
        Spk03[] = -Ecorr;
        Spk04[] = -Ecorr;

		Somega[]  = PO + PPO;
		Spomega[] = -EO;

	}
	//boundary({Sk01, Sk02, Sk03, Sk04,
	//Spk01, Spk02, Spk03, Spk04, Somega, Spomega});
}
#else
void compute_stress_omega_srcs(
	scalar Sk01, scalar Sk02, scalar Sk03,
	scalar Sk04, scalar Sk05, scalar Sk06,
	scalar Spk01, scalar Spk02, scalar Spk03,
	scalar Spk04, scalar Spk05, scalar Spk06,
scalar Somega, scalar Spomega
)
{
	tensor gradU[]; // grad(u)
	tensor S[], W[], tG[];
	vector df[], dk[], de[];

	t_rsm_centered_gradientv(u, gradU);
	t_rsm_centered_gradient(rhok, dk);
	t_rsm_centered_gradient(rhoe, de);

	foreach(){
		// compute dui_dxi
		double divU = gradU.x.x[] + gradU.y.y[] + gradU.z.z[];

		// Compute stress tensor S_{ij}
		double S_tr = (gradU.x.x[] + gradU.y.y[] + gradU.z.z[])/3.; // remove trace; divergence free
		S.x.x[] = gradU.x.x[];
		S.y.y[] = gradU.y.y[];
		S.z.z[] = gradU.z.z[];
		S.x.y[] = 0.5*(gradU.x.y[] + gradU.y.x[]);
		S.x.z[] = 0.5*(gradU.x.z[] + gradU.z.x[]);
		S.y.z[] = 0.5*(gradU.y.z[] + gradU.z.y[]);

		W.x.y[] = 0.5*(gradU.x.y[] - gradU.y.x[]);
		W.x.z[] = 0.5*(gradU.x.z[] - gradU.z.x[]);
		W.y.z[] = 0.5*(gradU.y.z[] - gradU.z.y[]);

		tG.x.x[] = S.x.x[] - S_tr;
		tG.y.y[] = S.y.y[] - S_tr;
		tG.z.z[] = S.z.z[] - S_tr;
		tG.x.y[] = S.x.y[];
		tG.x.z[] = S.x.z[];
		tG.y.z[] = S.y.z[];
		
		double alphaRho = 1.0;

#ifdef EMBED
		alphaRho *= clamp(cs[], 0.0, 1.0);
#endif

        // Compute R production and transportation
		double Pxx, Pxy, Pxz, Pyy, Pyz, Pzz;
		double Yxx, Yxy, Yxz, Yyy, Yyz, Yzz;
		double Dxx, Dxy, Dxz, Dyy, Dyz, Dzz;
        Pxx = (Rxx[]*gradU.x.x[] + Rxy[]*gradU.x.y[] + Rxz[]*gradU.x.z[])*2;
        Pxy =  Rxx[]*gradU.y.x[] + Rxy[]*gradU.y.y[] + Rxz[]*gradU.y.z[]
			 + Rxy[]*gradU.x.x[] + Ryy[]*gradU.x.y[] + Ryz[]*gradU.x.z[];
        Pxz =  Rxx[]*gradU.z.x[] + Rxy[]*gradU.z.y[] + Rxz[]*gradU.z.z[]
			 + Rxz[]*gradU.x.x[] + Ryz[]*gradU.x.y[] + Rzz[]*gradU.x.z[];
        Pyy = (Rxy[]*gradU.y.x[] + Ryy[]*gradU.y.y[] + Ryz[]*gradU.y.z[])*2;
        Pyz =  Rxy[]*gradU.z.x[] + Ryy[]*gradU.z.y[] + Ryz[]*gradU.z.z[]
			 + Rxz[]*gradU.y.x[] + Ryz[]*gradU.y.y[] + Rzz[]*gradU.y.z[];
        Pzz = (Rxz[]*gradU.z.x[] + Ryz[]*gradU.z.y[] + Rzz[]*gradU.z.z[])*2;

        Dxx = (Rxx[]*gradU.x.x[] + Rxy[]*gradU.y.x[] + Rxz[]*gradU.z.x[])*2;
        Dxy =  Rxx[]*gradU.x.y[] + Rxy[]*gradU.y.y[] + Rxz[]*gradU.z.y[]
			 + Rxy[]*gradU.x.x[] + Ryy[]*gradU.y.x[] + Ryz[]*gradU.z.x[];
        Dxz =  Rxx[]*gradU.x.z[] + Rxy[]*gradU.y.z[] + Rxz[]*gradU.z.z[]
			 + Rxz[]*gradU.x.x[] + Ryz[]*gradU.y.x[] + Rzz[]*gradU.z.x[];
        Dyy = (Rxy[]*gradU.x.y[] + Ryy[]*gradU.y.y[] + Ryz[]*gradU.z.y[])*2;
        Dyz =  Rxy[]*gradU.x.z[] + Ryy[]*gradU.y.z[] + Ryz[]*gradU.z.z[]
			 + Rxz[]*gradU.x.y[] + Ryz[]*gradU.y.y[] + Rzz[]*gradU.z.y[];
        Dzz = (Rxz[]*gradU.x.z[] + Ryz[]*gradU.y.z[] + Rzz[]*gradU.z.z[])*2;

        double P = 0.5*(Pxx + Pyy + Pzz);

        Yxx = kTurbConstants.beta_star*kTurbConstants.C1*rhoe[]*(2.0/3.0*rhok[]) - kTurbConstants.alpha_v*(Pxx - 2.0/3.0*P)
                    - kTurbConstants.beta_v*(Dxx - 2.0/3.0*P) - kTurbConstants.gamma_v*rhok[]*(tG.x.x[]);
        Yxy =-kTurbConstants.alpha_v*(Pxy) - kTurbConstants.beta_v*(Dxy) - kTurbConstants.gamma_v*rhok[]*(S.x.y[]);
        Yxz =-kTurbConstants.alpha_v*(Pxz) - kTurbConstants.beta_v*(Dxz) - kTurbConstants.gamma_v*rhok[]*(S.x.z[]);
        Yyy = kTurbConstants.beta_star*kTurbConstants.C1*rhoe[]*(2.0/3.0*rhok[]) - kTurbConstants.alpha_v*(Pyy - 2.0/3.0*P)
                    - kTurbConstants.beta_v*(Dyy - 2.0/3.0*P) - kTurbConstants.gamma_v*rhok[]*(tG.y.y[]);
        Yyz =-kTurbConstants.alpha_v*(Pyz) - kTurbConstants.beta_v*(Dyz) - kTurbConstants.gamma_v*rhok[]*(S.y.z[]);
        Yzz = kTurbConstants.beta_star*kTurbConstants.C1*rhoe[]*(2.0/3.0*rhok[]) - kTurbConstants.alpha_v*(Pzz - 2.0/3.0*P)
                    - kTurbConstants.beta_v*(Dzz - 2.0/3.0*P) - kTurbConstants.gamma_v*rhok[]*(tG.z.z[]);

#ifdef TWO_PHASE
		// Buoyancy correction, Pb = - rho p_b nu_t
		// buoyancy correction, p_b = \alpha_b^* G_i drhoidxj (drhodx = (rho2-rho1)*dfdx ) (issue: unbounded growth in air side?)
		double dfdx = (f[1, 0, 0]-f[-1, 0,  0])/2.0/Delta;
		double dfdy = (f[0, 1, 0]-f[0, -1,  0])/2.0/Delta;
		double dfdz = (f[0, 0, 1]-f[0,  0, -1])/2.0/Delta;

		double Pbcorr_xy = alphaRho * kTurbConstants.alpha_b_star * G.y * dfdx * nu_t[];
		double Pbcorr_yy = alphaRho * kTurbConstants.alpha_b_star * G.y * dfdy * nu_t[];
		double Pbcorr_yz = alphaRho * kTurbConstants.alpha_b_star * G.y * dfdz * nu_t[];
#else 
		double Pbcorr_xy = 0;
		double Pbcorr_yy = 0;
		double Pbcorr_yz = 0;
#endif

		double E = 2.0/3.0*kTurbConstants.beta_star*rhoe[]*rhok[];
		double Ecorr = kTurbConstants.beta_star*kTurbConstants.C1*rhoe[];

        // Compute the dissipation term
        // compute R_{ij}duidxj
		double RbyGradU = Rxx[] * gradU.x.x[] + Rxy[] * gradU.x.y[] + Rxz[] * gradU.x.z[] + \
				          Rxy[] * gradU.y.x[] + Ryy[] * gradU.y.y[] + Ryz[] * gradU.y.z[] + \
                          Rxz[] * gradU.z.x[] + Ryz[] * gradU.z.y[] + Rzz[] * gradU.z.z[] ;

        double PO = kTurbConstants.alpha*RbyGradU / (rhok[] + 1e-15) * rhoe[];
        double EO = kTurbConstants.beta*rhoe[];

		// Secondary generation of omega, PP_\omega = \sigma_d / \omega \frac{\partial_i}{k} \frac{\partial_i}{\omega}
		// compute dki_dxi domegai_dxi
		// notice the treatment to \nabla\omega/\omega (assume \omega is bounded to a positive value properly) to avoid a explicit division to \omega which could lead to numerical blow-up.
        double dkdw = dk.x[]*de.x[]+dk.y[]*de.y[]+dk.z[]*de.z[];
		double sigma_d = dkdw > 0.0? 0.125 : 0.0;
		double PPO = rhoe[] > kTurbConstants.omegaMin_ ? alphaRho * sigma_d * dkdw / rhoe[] : 0.0;

		// Update source terms
		Sk01[] = -Pxx-Yxx + E;
        Sk02[] = -Pxy-Yxy + Pbcorr_xy;
        Sk03[] = -Pxz-Yxz;
        Sk04[] = -Pyy-Yyy + E + Pbcorr_yy;
        Sk05[] = -Pyz-Yyz + Pbcorr_yz;
        Sk06[] = -Pzz-Yzz + E;

		Spk01[] = -Ecorr;
        Spk02[] = -Ecorr;
        Spk03[] = -Ecorr;
        Spk04[] = -Ecorr;
        Spk05[] = -Ecorr;
        Spk06[] = -Ecorr;

		Somega[]  = PO + PPO;
		Spomega[] = -EO;

	}
	//boundary({Sk01, Sk02, Sk03, Sk04, Sk05, Sk06,
	//Spk01, Spk02, Spk03, Spk04, Spk05, Spk06, Somega, Spomega});
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

#if dimension == 2
event turbulence_correction (i++) {
	if (kVerbose>1)
		fprintf(stdout,"[%d] Stress-omega: turbulence model correction.\n",i);

	// incorporate turbulence dissipation in k and omega equations
	face vector D_k[], D_omega[];
	compute_komega_mueff(D_k, D_omega);
	scalar Sk01[], Sk02[], Sk03[], Sk04[];
	scalar Spk01[], Spk02[], Spk03[], Spk04[];
    scalar src_omega[], sp_omega[];
	compute_stress_omega_srcs(Sk01, Sk02, Sk03, Sk04,
	Spk01, Spk02, Spk03, Spk04, src_omega, sp_omega);
	scalar rhov1[],rhov21[],rhov22[],rhov23[],rhov24[];

	// Source
	foreach()
	{
		for (scalar r in {rhov1, rhov21,rhov22,rhov23,rhov24}) {
			r[] = cm[]*alphaRho(f[]);
#if EMBED
			r[] *= max(cs[], 1e-14);
#endif
		}
        for (scalar s in {Sk01, Sk02, Sk03, Sk04,
			Spk01, Spk02, Spk03, Spk04, src_omega, sp_omega})
            s[] *= rhov1[];
	}
	//boundary({Rxx, Rxy, Ryy, Rzz, rhok, rhoe});

	// Advection
	if (kVerbose>1)
		fprintf(stdout,"[%d] Stress-omega: advect turbulence variables.\n",i);
	advection ((scalar *){Rxx, Rxy, Ryy, Rzz, rhoe}, uf, dt);
		
	// Diffusion
	if (kVerbose>1)
		fprintf(stdout,"[%d] Stress-Omega: diffuse turbulence variables.\n",i);

    bool mgFailed = 0;

	mgOmega = diffusion (rhoe, dt, D_omega, beta=sp_omega, r=src_omega, theta=rhov1);

	if (kVerbose>0)
	{
		fprintf (stderr,
				"MG convergence for rhoe reached after %d iterations, res: %g nrelax: %d\n",
				mgOmega.i, mgOmega.resa, mgOmega.nrelax);
	}

    scalar tmpR, tmpSU, tmpSP, tmpRHO;
	scalar *Rlist = (scalar *){Rxx, Rxy, Ryy, Rzz};
	scalar *SUlist = (scalar *){Sk01, Sk02, Sk03, Sk04};
	scalar *SPlist = (scalar *){Spk01, Spk02, Spk03, Spk04};
	scalar *RHOlist = (scalar *){rhov21,rhov22,rhov23,rhov24};

    for (tmpR, tmpSU, tmpSP, tmpRHO in Rlist, SUlist, SPlist, RHOlist)
    {
        mgstats mgR  = diffusion (tmpR, dt, D_k, beta=tmpSP , r=tmpSU , theta=tmpRHO);

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

    if( mgFailed || mgOmega.resa > 1.0)
    {
        debug_dump(i,(scalar *){f, nu_t, rhok, rhoe, Rxx, Rxy, Ryy, Rzz});
    }

	// Compute nu_t
	correct_nut();

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
#else
event turbulence_correction (i++) {
	if (kVerbose>1)
		fprintf(stdout,"[%d] Stress-omega: turbulence model correction.\n",i);

	// incorporate turbulence dissipation in k and omega equations
	face vector D_k[], D_omega[];
	compute_komega_mueff(D_k, D_omega);
	scalar Sk01[], Sk02[], Sk03[], Sk04[], Sk05[], Sk06[];
	scalar Spk01[], Spk02[], Spk03[], Spk04[], Spk05[], Spk06[];
    scalar src_omega[], sp_omega[];
	compute_stress_omega_srcs(Sk01, Sk02, Sk03, Sk04, Sk05, Sk06,
	Spk01, Spk02, Spk03, Spk04, Spk05, Spk06, src_omega, sp_omega);
	scalar rhov1[],rhov21[],rhov22[],rhov23[],rhov24[],rhov25[],rhov26[];

	// Source
	foreach()
	{
		for (scalar r in {rhov1, rhov21,rhov22,rhov23,rhov24,rhov25,rhov26}) {
			r[] = cm[]*alphaRho(f[]);
#if EMBED
			r[] *= max(cs[], 1e-14);
#endif
		}
        for (scalar s in {Sk01, Sk02, Sk03, Sk04, Sk05, Sk06,
	Spk01, Spk02, Spk03, Spk04, Spk05, Spk06, src_omega, sp_omega})
            s[] *= rhov1[];
	}
	//boundary({Rxx, Rxy, Rxz, Ryy, Ryz, Rzz, rhok, rhoe});

	// Advection
	if (kVerbose>1)
		fprintf(stdout,"[%d] Stress-omega: advect turbulence variables.\n",i);
	advection ((scalar *){Rxx, Rxy, Rxz, Ryy, Ryz, Rzz, rhoe}, uf, dt);
		
	// Diffusion
	if (kVerbose>1)
		fprintf(stdout,"[%d] Stress-Omega: diffuse turbulence variables.\n",i);

    bool mgFailed = 0;

	mgOmega = diffusion (rhoe, dt, D_omega, beta=sp_omega, r=src_omega, theta=rhov1);

	if (kVerbose>0)
	{
		fprintf (stderr,
				"MG convergence for rhoe reached after %d iterations, res: %g nrelax: %d\n",
				mgOmega.i, mgOmega.resa, mgOmega.nrelax);
	}

    scalar tmpR, tmpSU, tmpSP, tmpRHO;
	scalar *Rlist = (scalar *){Rxx, Rxy, Rxz, Ryy, Ryz, Rzz};
	scalar *SUlist = (scalar *){Sk01, Sk02, Sk03, Sk04, Sk05, Sk06};
	scalar *SPlist = (scalar *){Spk01, Spk02, Spk03, Spk04, Spk05, Spk06};
	scalar *RHOlist = (scalar *){rhov21,rhov22,rhov23,rhov24,rhov25,rhov26};

    for (tmpR, tmpSU, tmpSP, tmpRHO in Rlist, SUlist, SPlist, RHOlist)
    {
        mgstats mgR  = diffusion (tmpR, dt, D_k, beta=tmpSP , r=tmpSU , theta=tmpRHO);

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

    if( mgFailed || mgOmega.resa > 1.0)
    {
        debug_dump(i,(scalar *){f, nu_t, rhok, rhoe, Rxx, Rxy, Rxz, Ryy, Ryz, Rzz});
    }

	// Compute nu_t
	correct_nut();

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
#endif

#if TREE
// mesh adaptation happens before solving viscous term
event adapt (i++) {
	if (kVerbose>1)
		fprintf(stdout,"[%d] Stress-omega: Correct nu_t after regridding.\n",i);
	correct_nut();
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

