#ifndef RSM_COMMON_H
#define RSM_COMMON_H
#include "levelset_common.h"

face vector mu_t[]; // turbulent viscosity
scalar nu_t[]; // Cell Centered diffusivity

#ifndef TWO_PHASE
scalar f[]; // dummy vof field
#ifndef RHO
#define RHO (1.0)
#endif
#ifndef MU
#define MU (1.0e-6)
#endif
#define alphaRho(f) (RHO)
#define alphaMu(f) (MU)
#define alphaMuEff(f) (MU)
#define alphaNuEff(f) (MU/RHO)
#else
#define alphaMu(f) (mu(f))
#define alphaRho(f) (rho(f))
#define alphaMuEff(f) (mu(f))
#define alphaNuEff(f) (mu(f)/rho(f))
#endif

// Abstract interface to implement
void correct_nut();
void set_rans_init();
void rsm_model_debug_dump(int i, scalar* plist);
void rsm_model_output(scalar**);
void rsm_model_output_compact(scalar**);
void refresh_rsm_model_by_json(struct rsm_JSONFile_t*);

// Utility functions
// note: u is weighted by fm
void t_rsm_centered_gradient (scalar p, vector g)
{

  /**
  We first compute a face field $\mathbf{g}_f$. */

  face vector gf[];
  foreach_face()
    gf.x[] = fm.x[]*(p[] - p[-1])/Delta;

  /**
  We average these face values to obtain the centered gradient field. */

  trash ({g});
  foreach()
    foreach_dimension()
      g.x[] = (gf.x[] + gf.x[1])/(fm.x[] + fm.x[1] + SEPS);
}

void t_rsm_centered_gradientv (vector p, tensor g){
	foreach_dimension()
		t_rsm_centered_gradient (p.x, g.x);
}

void t_rsm_biased_gradient (scalar p, vector g, bool inverse)
{
  /**
  We first compute a face field $\mathbf{g}_f$. */

  face vector gf[];
  foreach_face()
    gf.x[] = fm.x[]*(p[] - p[-1])/Delta;

  /**
  We average these face values to obtain the centered gradient field. */

  trash ({g});
  foreach()
    foreach_dimension() {
		g.x[] = (gf.x[]*(inverse ? (1.-f[-1]):f[-1]) + gf.x[1]*(inverse ? (1.-f[1]):f[1]))
		/(fm.x[] + fm.x[1] + SEPS)
		/( (inverse ? (1.-f[-1]):f[-1]) + (inverse ? (1.-f[1]):f[1]) + 1e-15);
	}
}

void t_rsm_biased_gradientv (vector p, tensor g, bool inverse){
	foreach_dimension()
		t_rsm_biased_gradient (p.x, g.x, inverse);
}

void bound(scalar p, double pMin)
{
	scalar q[];

	foreach(){
		p[] = p[] > pMin ? p[] : pMin;
		q[] = p[];
	}

	//boundary({q});

	foreach(){
		if(q[] == pMin)
		{
			int nb = 0;
			double nbVals = 0.0;
			foreach_neighbor(1)
			{
				if (q[] > pMin)
				{
					nb++;
					nbVals += q[];
				}
			}
			if(nb > 0)
				p[] = nbVals / nb;
		}
	}
}

void bound2(scalar p, double pMin, double pMax)
{
	foreach(){
		p[] = p[] > pMin ? p[] : pMin;
		p[] = p[] < pMax ? p[] : pMax;
	}
	foreach(){
		if(p[] == pMin || p[] == pMax)
		{
			int nb = 0;
			double nbVals = 0.0;
			foreach_neighbor(1)
			{
				if (p[] > pMin && p[] < pMax)
				{
					nb++;
					nbVals += p[];
				}
			}
			if(nb > 0)
				p[] = nbVals / nb;
		}
	}
}

// Include implementation of rans model
#ifdef RSM_KOMEGA
	#include "komega.h"
	#define RSM_IMPL
#endif

#ifdef RSM_STRESS_OMEGA
	#include "stress_omega.h"
	#define RSM_IMPL
#endif

// Define dummy functions
#ifndef RSM_IMPL
	void correct_nut(){}
	void set_rans_init(){}
	void rsm_model_debug_dump(int i, scalar* plist){}
	void rsm_model_output(scalar** plist){}
	void rsm_model_output_compact(scalar** plist){}
	void refresh_rsm_model_by_json(struct rsm_JSONFile_t* jf){}
#endif

#endif

