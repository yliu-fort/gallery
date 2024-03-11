# config.json requires following fields
"verbose":0,
"size": [3.79, 83.79, -0.8, 79.2],
"waterlevel": 0, // No implementation yet. Do not recommend to modify this.
"g": 9.81,
"rho1": 1025,
"rho2": 1.225,
"mu1": 0.00089,
"mu2": 1.74E-5,
"refine_init": 7,
"adapt_refine_max": 11,
"adapt_refine_min": 5,
"adapt_utol": 0.002,

"t0": 0.0,
"tmax": 130.0,
"ramptime": 1.0, // set to 0 if generating a spatial initial condition.
"cfl": 0.1,
"drift_velocity": 0,

# if RSM model is enabled then
"k0": 1E-4,
"omega0": 421.0

# wave_profile.json requires following fields
{
  "Hs": 0.01586376658246301,
  "w0": 5.723745802214097,
  "k0": 3.3395785941247196,
  "bsk": {
    "nks1": 1,
    "nks2": 0,
    "nks3": 0,
    "nks": 1,
    "wavePeriods": [7.1651308522830375],
    "waveLengths": [80.156250000000711],
    "waveHeights": [0],
    "wavePhases": [-21.812217565446748],
    "waveDirs": [0],
    "t_start": 0.0,
    "t_end": 20.0
  },

  "globals": {
    
	
	# Optional fields to determine the size, t0, tmax and depth
	# These fields have higher priority against config.json
	"use_default_grid": true, // set to false to use x_min & x_max below.
    "x_min": -40,
    "x_max": 40,
    "nx": 513,
    "y_min": -1,
    "y_max": 1,
    "ny": 3,
    "nt": 200,
	"duration": [0, 20],
	"depth": 0.8,
  }
}

# wave_gauge.json requires following fields
{
	"type" : "line_distributed",
	"x0" : 0.0,
	"x1" : 1.0,
  "y0" : 0.6, // (todo) implement y1 and z1
	"z0" : 0.0,
	"absolute_pos" : false,
	"num_probes" : 513
}
or
{
	"type" : "user_defined",
	"pos_wg_x" : [3.79, 6.82, 9.82],
  "y0" : 0.6, // (todo) implement pos_wg_y and pos_wg_z, now just uniform
	"z0" : 0.0,
	"absolute_pos" : true,
	"num_probes" : 3
}

# Default parameters that is determined during compilation
static const double SURFACE_OUTPUT_PERIOD = 0.01; // second
static const double OUTPUT_PERIOD=0.1; // second
static const double AUTOSAVE_PERIOD=1.0; // second
static const double CHECKPOINT_PERIOD=10.0; // second
static const double MAXIMUM_RUNTIME = 300.0; // second
static const double SUBGRID_DROPLET_REMOVAL_PERIOD=1.0; // second

# Compile flags
-DRSM_* : RSM model related flags
-DRSM_KOMEGA : turn on k-omega model
-DRSM_STRESS_OMEGA : turn on stress-omega model
-DRSM_NO_INTERFACE_DIFFUSION : turn off turbulence viscosity near interface
-DNO_INLET : set to periodic domain / non-periodic domain
-DSTOKES : set to analytic stokes wave / user defined wave profile