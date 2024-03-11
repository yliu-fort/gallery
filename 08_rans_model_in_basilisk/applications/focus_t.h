#ifndef STOKES
double wave (double x, double y, double t); // elevation
double u_x (double x, double y, double t); //velocity_x
double u_y (double x, double y, double t); //velocity_y
double u_z (double x, double y, double t); //velocity_z

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //
double etaMulti
(
    const double H,
    const double Kx,
    const double xi,
    const double Ky,
    const double yi,
    const double omega,
    const double ti,
    const double phase
)
{
    double phaseTot = Kx*xi + Ky*yi - omega*ti + phase;
    return H*cos(phaseTot);
}

double uMultiDirec_x
(
    const double irregH,
    const double waveOmega,
    const double pha,
    const double waveKs,
    const double zz,
    const double hh,
    const double irregDir
)
{
    const double ksh = waveKs*hh;
    const double ksz = waveKs*zz;

    return  irregH*waveOmega*cos(pha)*(coshl(ksz)/sinhl(ksh))*cos(irregDir);
}

double uMultiDirec_z
(
    const double irregH,
    const double waveOmega,
    const double pha,
    const double waveKs,
    const double zz,
    const double hh,
    const double irregDir
)
{
    const double ksh = waveKs*hh;
    const double ksz = waveKs*zz;

    return  irregH*waveOmega*cos(pha)*(coshl(ksz)/sinhl(ksh))*sin(irregDir);
}

double uMultiDirec_y
(
    const double irregH,
    const double waveOmega,
    const double pha,
    const double waveKs,
    const double zz,
    const double hh,
    const double irregDir
)
{
    const double ksh = waveKs*hh;
    const double ksz = waveKs*zz;
	
    return irregH*waveOmega*sin(pha)*(sinhl(ksz)/sinhl(ksh));
}

double wave (double x, double y, double t)
{
	if (t > kMultiwaveInlet.t_end_)
	{
		return 0.0;
	}

    double eta = 0;

    for ( int ii = 0; ii < kMultiwaveInlet.nks_; ii++)
    {

        double waveOrder = ii < kMultiwaveInlet.nks1_? 1 : (ii < (kMultiwaveInlet.nks1_ + kMultiwaveInlet.nks2_)? 2 : 3);
        double waveKs = 2*M_PI/kMultiwaveInlet.waveLengths_[ii];
        double waveOmegas = 2*M_PI/kMultiwaveInlet.wavePeriods_[ii];

            eta += etaMulti
                (
                    kMultiwaveInlet.waveHeights_[ii],
                    waveKs*cos(kMultiwaveInlet.waveDirs_[ii]),
                    x,
                    waveKs*sin(kMultiwaveInlet.waveDirs_[ii]),
                    0,
                    waveOmegas,
                    t,
                    kMultiwaveInlet.wavePhases_[ii] + kWaveData.phase_shift * waveOrder
                );
    }

    return eta;
}

double u_x (double x, double y, double t)
{
    if (t > kMultiwaveInlet.t_end_)
	{
		return 0.0;
	}

    if (y > kWaveData.waterdepth)
    {
        return 0.0;
    }

    double u = 0.0;

    for ( int ii = 0; ii < kMultiwaveInlet.nks1_; ii++)
    {
        double waveKs = 2.0*M_PI/kMultiwaveInlet.waveLengths_[ii];
        double waveOmegas = 2.0*M_PI/kMultiwaveInlet.wavePeriods_[ii];

        double phaseTot =
            waveKs*x*cos(kMultiwaveInlet.waveDirs_[ii])
            //+ waveKs*z*sin(kMultiwaveInlet.waveDirs_[ii])
            - waveOmegas*t
            + kMultiwaveInlet.wavePhases_[ii]
            + kWaveData.phase_shift;
            

        const double Uf = uMultiDirec_x
        (
            kMultiwaveInlet.waveHeights_[ii],
            waveOmegas,
            phaseTot,
            waveKs,
            y+kWaveData.waterdepth,
            kWaveData.waterdepth,
            kMultiwaveInlet.waveDirs_[ii]
        );

        u += Uf;
    }

    return u;
}

double u_y (double x, double y, double t)
{
    if ( t > kMultiwaveInlet.t_end_)
	{
		return 0.0;
	}

    if (y > kWaveData.waterdepth)
    {
        return 0.0;
    }

    double u = 0.0;

    for ( int ii = 0; ii < kMultiwaveInlet.nks1_; ii++)
    {
        double waveKs = 2.0*M_PI/kMultiwaveInlet.waveLengths_[ii];
        double waveOmegas = 2.0*M_PI/kMultiwaveInlet.wavePeriods_[ii];

        double phaseTot =
            waveKs*x*cos(kMultiwaveInlet.waveDirs_[ii])
            //+ waveKs*z*sin(kWaveData.waveDirs_[ii])
            - waveOmegas*t
            + kMultiwaveInlet.wavePhases_[ii]
            + kWaveData.phase_shift;
            

        const double Uf = uMultiDirec_y
        (
            kMultiwaveInlet.waveHeights_[ii],
            waveOmegas,
            phaseTot,
            waveKs,
            y+kWaveData.waterdepth,
            kWaveData.waterdepth,
            kMultiwaveInlet.waveDirs_[ii]
        );

        u += Uf;
    }

    return u;
}

double u_z (double x, double y, double t)
{
    if (t > kMultiwaveInlet.t_end_)
	{
		return 0.0;
	}

    double u = 0.0;

    for ( int ii = 0; ii < kMultiwaveInlet.nks1_; ii++)
    {
        double waveKs = 2.0*M_PI/kMultiwaveInlet.waveLengths_[ii];
        double waveOmegas = 2.0*M_PI/kMultiwaveInlet.wavePeriods_[ii];

        double phaseTot =
            waveKs*x*cos(kMultiwaveInlet.waveDirs_[ii])
            //+ waveKs*z*sin(kMultiwaveInlet.waveDirs_[ii])
            - waveOmegas*t
            + kMultiwaveInlet.wavePhases_[ii]
            + kWaveData.phase_shift;
            

        const double Uf = uMultiDirec_z
        (
            kMultiwaveInlet.waveHeights_[ii],
            waveOmegas,
            phaseTot,
            waveKs,
            y+kWaveData.waterdepth,
            kWaveData.waterdepth,
            kMultiwaveInlet.waveDirs_[ii]
        );

        u += Uf;
    }

    return u;
}
#endif
