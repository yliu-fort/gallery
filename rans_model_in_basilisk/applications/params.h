#ifndef WAVE_PARAMS_H
#define WAVE_PARAMS_H

#include "ray_math.h"

struct WaveProperties kWaveData;
struct MultidirectionInletProfile kMultiwaveInlet;
struct WaveGaugeProfile kWaveGauge;

struct WaveGaugeProfile {
    Ray* rays;
    size_t num_probes;
};

struct MultidirectionInletProfile {
    double t_start_;
    double t_end_;
    int nks_;
    int nks1_;
    int nks2_;
    int nks3_;
    double *wavePeriods_;
    double *waveLengths_;
    double *waveHeights_;
    double *wavePhases_;
    double *waveDirs_;
};

struct WaveProperties {
	// Wave configuration
	double t0; // simulation start time, s
	double tmax; // simulation end time, s
    double dtmax; // simulation max time interval, s
	double *size; // waveflume dimensions [x_min, x_max, y_min, y_max, z_min, z_max], m
  	double waterlevel; // vertical shift of reference water level relative to domain definition (default 0)
  	double g; // gravity, m/s^2 (default 9.81)
  	double rho1; // density of water, kg/m^3
  	double rho2; // density of air, kg/m^3
	double mu1; // viscosity of water, kg/m^3
	double mu2; // viscosity of air, kg/m^3
    double re; // Reynolds number water
    double mu_ratio; // ratio of viscosity of water/viscosity of air
	int refine_init; // refinment level at t = 0 (default 7)
	int adapt_refine_max; // maximum level of refinement (default 14)
	int adapt_refine_min; // minimum level of refinement (default 5)
	double adapt_utol; // refinement criterion: velocity tolerance (default 0.005)
	double Hs; // Significant wave height, m
	double T0; // wave period, s
	double L0; // wave number, rad/m
	double waterdepth; // waterdepth, m
	double ramptime; // ramptime, s (default 1.0s)
	double cfl; // maximum allowed courant number (default 0.1)
    double cfl_ramptime; // courant number ramp-up time (default 2.0)
	double drift_velocity; // drift velocity, m/s (default is 0)
    double phase_shift;
};

void init_wavedata_default() {
    kWaveData.t0 = 0.0;
    kWaveData.tmax = 100.0;
    kWaveData.dtmax = 0.005;
    if(!kWaveData.size) {kWaveData.size = calloc(6, sizeof(double));}
#ifdef STOKES
    kWaveData.size[0] = 0.0;
    kWaveData.size[1] = 12.4914;
    kWaveData.size[2] = -7.33;
    kWaveData.size[3] = 5.1614;
    kWaveData.size[4] = -1.0;
    kWaveData.size[5] = 1.0;
#else
    kWaveData.size[0] = 3.790000;
    kWaveData.size[1] = 83.790000;
    kWaveData.size[2] = -0.800000;
    kWaveData.size[3] = 1.000000;
    kWaveData.size[4] = -1.0;
    kWaveData.size[5] = 1.0;
#endif
    kWaveData.waterlevel = 0;
    kWaveData.g = 9.81;
    kWaveData.rho1 = 1025.0;
    kWaveData.rho2 = 1.225;
    kWaveData.mu1 = 8.9e-4;
    kWaveData.mu2 = 17.4e-6;
    kWaveData.re = 1000000;
    kWaveData.mu_ratio = 50;
#if dimension == 2
    kWaveData.refine_init = 7;
    kWaveData.adapt_refine_max = 14;
    kWaveData.adapt_refine_min = 5;
#elif dimension == 3
    kWaveData.refine_init = 5;
    kWaveData.adapt_refine_max = 7;
    kWaveData.adapt_refine_min = 3;
#endif
    kWaveData.adapt_utol = 0.002;
#ifdef STOKES
    kWaveData.Hs = 0.2;
    kWaveData.T0 = 2.0;
    kWaveData.L0 = 6.28;
    kWaveData.waterdepth = 7.33;
    kWaveData.ramptime = 0;
#else
    kWaveData.Hs = 0.02;
    kWaveData.T0 = 1.0;
    kWaveData.L0 = 1.0;
    kWaveData.waterdepth = 0.8;
    kWaveData.ramptime = 1.0;
#endif
    kWaveData.cfl = 0.1;
    kWaveData.cfl_ramptime = 2.0; // courant number ramp-up time (default 2.0)
    kWaveData.drift_velocity = 0.0;

    kWaveData.phase_shift = 0.0;
}

void PrintWaveProperties(struct WaveProperties *args)
{
	  //printf ("NSUB = %f\n", args->nsub);
	  //printf ("dataset_path = %s\n",args->dataset_path);
	  //printf ("dataset_name = %s\n",args->dataset_name);

	  // Wave configuration
	  printf ("time_start = %f\n",args->t0);
	  printf ("time_end = %f\n",args->tmax);
      printf ("max_time_interval = %f\n",args->dtmax);
#if dimension == 2
	  printf ("domain_size = (%f, %f, %f, %f)\n",args->size[0],args->size[1],args->size[2],args->size[3]);
#elif dimension == 3
	  printf ("domain_size = (%f, %f, %f, %f, %f, %f)\n",
			  args->size[0],args->size[1],args->size[2],args->size[3],args->size[4],args->size[5]);
#endif
	  printf ("ref_level = %f\n",args->waterlevel);
	  printf ("gravity = %f\n",args->g);
	  printf ("rho_water = %f\n",args->rho1);
	  printf ("rho_air = %f\n",args->rho2);
	  printf ("mu_water = %f\n",args->mu1);
	  printf ("mu_air = %f\n",args->mu2);
      printf ("re = %f\n",args->re);
      printf ("mu_ratio = %f\n",args->mu_ratio);
	  printf ("initial_refinement_level = %d\n",args->refine_init);
	  printf ("max_refinement_level = %d\n",args->adapt_refine_max);
	  printf ("min_refinement_level = %d\n",args->adapt_refine_min);
	  printf ("adaptive_tolerence_u = %f\n",args->adapt_utol);
	  printf ("significant_wave_height = %f\n",args->Hs);
	  printf ("characteristic_period = %f\n",args->T0);
	  printf ("characteristic_length = %f\n",args->L0);
	  printf ("water_depth = %f\n",args->waterdepth);
	  printf ("ramp_time = %f\n",args->ramptime);
	  printf ("courant_number = %f\n",args->cfl);
      printf ("CFL ramp-up time = %f\n",args->cfl_ramptime);
	  printf ("drift_velocity = %f\n",args->drift_velocity);
      printf ("phase_shift = %f\n",args->phase_shift);
}

// *** //
void refresh_waveinput_by_json(struct rsm_JSONFile_t *jf, struct rsm_JSONFile_t *jf2)
{
    const cJSON *system_config = jf->json;
    rsm_print_json(system_config);

    read_double(&(kWaveData.waterlevel), "waterlevel", system_config);
    read_double(&(kWaveData.g), "g", system_config);
    read_double(&(kWaveData.rho1), "rho1", system_config);
    read_double(&(kWaveData.rho2), "rho2", system_config);
    read_double(&(kWaveData.mu1), "mu1", system_config);
    read_double(&(kWaveData.mu2), "mu2", system_config);
    read_double(&(kWaveData.re), "re", system_config);
    read_double(&(kWaveData.mu_ratio), "mu_ratio", system_config);
    read_int(&(kWaveData.refine_init), "refine_init", system_config);
    read_int(&(kWaveData.adapt_refine_max), "adapt_refine_max", system_config);
    read_int(&(kWaveData.adapt_refine_min), "adapt_refine_min", system_config);
    read_double(&(kWaveData.adapt_utol), "adapt_utol", system_config);
    read_double(&(kWaveData.ramptime), "ramptime", system_config);
    read_double(&(kWaveData.cfl), "cfl", system_config);
    read_double(&(kWaveData.cfl_ramptime), "cfl_ramptime", system_config);
    read_double(&(kWaveData.drift_velocity), "drift_velocity", system_config);
    read_double(&(kWaveData.phase_shift), "phase_shift", system_config);
    read_double(&(kWaveData.t0), "t0", system_config);
    read_double(&(kWaveData.tmax), "tmax", system_config);
    read_double(&(kWaveData.dtmax), "dtmax", system_config);
	
	const cJSON *wave_config = jf2->json;
    rsm_print_json(wave_config);
	const cJSON *wave_global_config = cJSON_GetObjectItemCaseSensitive(wave_config, "globals");
    if(wave_global_config){rsm_print_json(wave_global_config);}

	read_double(&(kWaveData.Hs), "Hs", wave_config);
    read_double(&(kWaveData.T0), "w0", wave_config);
	kWaveData.T0 = kWaveData.T0 > 0 ? 2.0 * 3.1415926 / kWaveData.T0 : 0.0;
    read_double(&(kWaveData.L0), "k0", wave_config);
	kWaveData.L0 = kWaveData.L0 > 0 ? 2.0 * 3.1415926 / kWaveData.L0 : 0.0;

    bool adim_input = false;
    read_bool(&adim_input, "adim_input", wave_global_config);
    double L0 = 1;
    double T0 = 1;
    double D0 = 1;
    if(adim_input)
    {
        double eps0 = 0.35;
        double k0 = 2.0 * 3.1415926 / kWaveData.L0;
        double w0 = 2.0 * 3.1415926 / kWaveData.T0;
        L0 = 1.0/(2.0*sq(eps0)*k0);
        T0 = 1.0/(w0*eps0);
        D0 = k0;
    }

	double *duration = NULL;
	read_num_array_1d(&duration, "duration", wave_global_config);
	if(duration != NULL)
	{
		kWaveData.tmax = duration[1]*T0;
        kWaveData.t0 += duration[0]*T0;
        free(duration);
        duration = NULL;
	}

    bool use_default_grid = true;
    read_bool(&use_default_grid, "use_default_grid", wave_global_config);
    if(use_default_grid)
    {
        read_num_array_1d(&(kWaveData.size), "size", system_config);
        read_double(&(kWaveData.waterdepth), "waterdepth", system_config);
    }else{
        read_double(&(kWaveData.waterdepth), "depth", wave_global_config);
        kWaveData.waterdepth /= D0;
        if(!kWaveData.size)
            kWaveData.size = calloc(6, sizeof(double));
        read_double(&(kWaveData.size[0]), "x_min", wave_global_config);
        read_double(&(kWaveData.size[1]), "x_max", wave_global_config);
#if dimension == 3
        read_double(&(kWaveData.size[4]), "y_min", wave_global_config);
        read_double(&(kWaveData.size[5]), "y_max", wave_global_config);
#else
        kWaveData.size[4] = -1;
        kWaveData.size[5] = 1;
#endif
        kWaveData.size[0] *= L0*(513.0/512.0); //(!!!WARNING!!!TODO) temorary fix to discontinuous IC
        kWaveData.size[1] *= L0*(513.0/512.0); //(!!!WARNING!!!TODO) temorary fix to discontinuous IC
        kWaveData.size[2] = -kWaveData.waterdepth;
        kWaveData.size[3] = kWaveData.size[1] - kWaveData.size[0] - kWaveData.waterdepth;
        kWaveData.size[4] *= L0*(513.0/512.0); //(!!!WARNING!!!TODO) temorary fix to discontinuous IC
        kWaveData.size[5] *= L0*(513.0/512.0); //(!!!WARNING!!!TODO) temorary fix to discontinuous IC
    }

    bool use_reynolds = false;
    read_bool(&(use_reynolds), "use_reynolds", system_config);
    if(use_reynolds)
    {
        kWaveData.mu1 = sqrt(kWaveData.g*cube(kWaveData.L0))/kWaveData.re*kWaveData.rho1;
        kWaveData.mu2 = kWaveData.mu1/kWaveData.mu_ratio;
    }
}

void refresh_multiwave_inlet_by_json(struct rsm_JSONFile_t *jf)
{
    const cJSON *wave_config = jf->json;
    const cJSON *multiinlet_config = cJSON_GetObjectItemCaseSensitive(wave_config, "bsk");
    if(!multiinlet_config) 
    {
        // search in root level
        multiinlet_config = wave_config;
    }
    rsm_print_json(multiinlet_config);
    read_int(&(kMultiwaveInlet.nks_), "nks", multiinlet_config);
    read_int(&(kMultiwaveInlet.nks1_), "nks1", multiinlet_config);
    read_int(&(kMultiwaveInlet.nks2_), "nks2", multiinlet_config);
    read_int(&(kMultiwaveInlet.nks3_), "nks3", multiinlet_config);
    read_double(&(kMultiwaveInlet.t_start_), "t_start", multiinlet_config);
    read_double(&(kMultiwaveInlet.t_end_), "t_end", multiinlet_config);
    read_num_array_1d(&(kMultiwaveInlet.wavePeriods_), "wavePeriods", multiinlet_config);
    read_num_array_1d(&(kMultiwaveInlet.waveLengths_), "waveLengths", multiinlet_config);
    read_num_array_1d(&(kMultiwaveInlet.waveHeights_), "waveHeights", multiinlet_config);
    read_num_array_1d(&(kMultiwaveInlet.wavePhases_), "wavePhases", multiinlet_config);
    read_num_array_1d(&(kMultiwaveInlet.waveDirs_), "waveDirs", multiinlet_config);

    //for(int i = 0; i < kMultiwaveInlet.nks_; i++)
    //    kMultiwaveInlet.wavePhases_[i] *= -1;
}

void refresh_wave_gauge_by_json(struct rsm_JSONFile_t *jf){
    const cJSON *wave_gauge_config = jf->json;
    if(!wave_gauge_config) {return;}

    char probe_type[255];
    read_char(probe_type, "type", wave_gauge_config);
    if(strcmp(probe_type, "line_distributed") == 0)
    {
            // statements
            int num_probes;
            double x0, x1, y0, z0;
            bool absolute_pos=false;

            read_double(&x0, "x0", wave_gauge_config);
            read_double(&x1, "x1", wave_gauge_config);
            read_double(&y0, "y0", wave_gauge_config);
            read_double(&z0, "z0", wave_gauge_config);
            read_int(&num_probes, "num_probes", wave_gauge_config);
            read_bool(&absolute_pos, "absolute_pos", wave_gauge_config);

            kWaveGauge.num_probes = num_probes;
            if(kWaveGauge.rays){free(kWaveGauge.rays);}
            kWaveGauge.rays = calloc(kWaveGauge.num_probes, sizeof(Ray));

            double xmin = (absolute_pos ? 0 : kWaveData.size[0]) + x0 * (absolute_pos ? 1 : (kWaveData.size[1] - kWaveData.size[0]));
            double dx = (x1 - x0) / (kWaveGauge.num_probes - 1) * (absolute_pos ? 1 : (kWaveData.size[1] - kWaveData.size[0]));
            for(int i = 0; i < kWaveGauge.num_probes; i++)
            {
                kWaveGauge.rays[i].orig.x = xmin + dx*i;
                kWaveGauge.rays[i].orig.y = y0;
                kWaveGauge.rays[i].orig.z = z0;
                kWaveGauge.rays[i].direction.x = 0.0;
                kWaveGauge.rays[i].direction.y =-1.0;
                kWaveGauge.rays[i].direction.z = 0.0;
            }
    }else if(strcmp(probe_type, "xz_plane_distributed") == 0){
            // statements
            double* num_probes = NULL; // nx, nz
            double x0, x1, y0, z0, z1;
            bool absolute_pos=false;

            read_double(&x0, "x0", wave_gauge_config);
            read_double(&x1, "x1", wave_gauge_config);
            read_double(&y0, "y0", wave_gauge_config);
            read_double(&z0, "z0", wave_gauge_config);
            read_double(&z1, "z1", wave_gauge_config);
            read_num_array_1d(&num_probes, "num_probes", wave_gauge_config);
            read_bool(&absolute_pos, "absolute_pos", wave_gauge_config);

            kWaveGauge.num_probes = (int)(num_probes[0])*(int)(num_probes[1]);
            if(kWaveGauge.rays){free(kWaveGauge.rays);}
            kWaveGauge.rays = calloc(kWaveGauge.num_probes, sizeof(Ray));

            double xmin = (absolute_pos ? 0 : kWaveData.size[0]) + x0 * (absolute_pos ? 1 : (kWaveData.size[1] - kWaveData.size[0]));
            double dx = (x1 - x0) / (num_probes[0] - 1) * (absolute_pos ? 1 : (kWaveData.size[1] - kWaveData.size[0]));
            double zmin = (absolute_pos ? 0 : kWaveData.size[4]) + z0 * (absolute_pos ? 1 : (kWaveData.size[5] - kWaveData.size[4]));
            double dz = (z1 - z0) / (num_probes[1] - 1) * (absolute_pos ? 1 : (kWaveData.size[5] - kWaveData.size[4]));

            int k = 0;
            for(int i = 0; i < (int)(num_probes[0]); i++)
            {
                for(int j = 0; j < (int)(num_probes[1]); j++)
                {
                    kWaveGauge.rays[k].orig.x = xmin + dx*i;
                    kWaveGauge.rays[k].orig.y = y0;
                    kWaveGauge.rays[k].orig.z = zmin + dz*j;
                    kWaveGauge.rays[k].direction.x = 0.0;
                    kWaveGauge.rays[k].direction.y =-1.0;
                    kWaveGauge.rays[k].direction.z = 0.0;
                    k++;
                }
            }
    }else if(strcmp(probe_type, "user_defined") == 0){
            // statements
            int num_probes;
            double y0, z0;
            double* pos_wg_x = NULL;
            bool absolute_pos=false;

            read_double(&y0, "y0", wave_gauge_config);
            read_double(&z0, "z0", wave_gauge_config);
            read_num_array_1d(&pos_wg_x, "pos_wg_x", wave_gauge_config);
            read_int(&num_probes, "num_probes", wave_gauge_config);
            read_bool(&absolute_pos, "absolute_pos", wave_gauge_config);

            kWaveGauge.num_probes = num_probes;
            if(kWaveGauge.rays){free(kWaveGauge.rays);}
            kWaveGauge.rays = calloc(kWaveGauge.num_probes, sizeof(Ray));

            for(int i = 0; i < kWaveGauge.num_probes; i++)
            {
                kWaveGauge.rays[i].orig.x = (absolute_pos ? 0 : kWaveData.size[0]) + pos_wg_x[i] * (absolute_pos ? 1 : (kWaveData.size[1] - kWaveData.size[0]));
                kWaveGauge.rays[i].orig.y = y0;
                kWaveGauge.rays[i].orig.z = z0;
                kWaveGauge.rays[i].direction.x = 0.0;
                kWaveGauge.rays[i].direction.y =-1.0;
                kWaveGauge.rays[i].direction.z = 0.0;
            }

            if(pos_wg_x){free(pos_wg_x);}
    }else{
        // default statements
        fprintf(stderr, "Unknown wave gauge type!"); 
    }
}


#endif
