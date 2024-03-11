/**
# Third-order Stokes wave

These functions return the shape of a third-order Stokes wave with the
wavenumber and steepness given by the parameters $ak$ and $k_$. */

#ifdef STOKES
 double wave (double x, double y, double t)
 {
     double ak_ = kWaveData.Hs/2.0*(2.0*M_PI/kWaveData.L0);
     double k_ = 2.0*M_PI/kWaveData.L0;
     double g_ = kWaveData.g;

     double H = ak_/k_*2.0;
     double h_ = kWaveData.waterdepth;
     double omega = sqrt(g_*k_*tanh(k_*h_));

     double phaseTot = omega*t - k_*x;
  
     return         H/2.0*cos(phaseTot)
                   +1.0/16.0*k_*sq(H)
                   *(
                        3.0/cube(tanh(k_*h_))
                      - 1.0/tanh(k_*h_)
                    )*cos(2.0*phaseTot);
 }
  
 double u_x (double x, double y, double t)
 {
    if (y > kWaveData.waterdepth)
    {
        return 0.0;
    }
     double ak_ = kWaveData.Hs/2.0*(2.0*M_PI/kWaveData.L0);
     double k_ = 2.0*M_PI/kWaveData.L0;
     double g_ = kWaveData.g;

     double H = ak_/k_*2.0;
     double h = kWaveData.waterdepth;
     double omega = sqrt(g_*k_*tanh(k_*h));

     double cel = omega/k_;

     double phaseTot = - k_*x + omega*t;

     // First order contribution
     double Uhorz = H/2.0*omega*
                   cosh(k_*(y + h))/sinh(k_*h) *
                   cos(phaseTot);

    // Second order contribution
    Uhorz += 3.0/16.0*cel*sq(k_*H)*cosh(2.0*k_*(y + h))
            /pow(sinh(k_*h),4.0)*cos(2.0*phaseTot)
             - 1.0/8.0*abs(g_)*sq(H)/(cel*h);
     return Uhorz;
 }
  
 double u_y (double x, double y, double t)
 {
    if (y > kWaveData.waterdepth)
    {
        return 0.0;
    }
     double ak_ = kWaveData.Hs/2.0*(2.0*M_PI/kWaveData.L0);
     double k_ = 2.0*M_PI/kWaveData.L0;
     double g_ = kWaveData.g;

     double H = ak_/k_*2.0;
     double h = kWaveData.waterdepth;
     double omega = sqrt(g_*k_*tanh(k_*h));

     double cel = omega/k_;

     double phaseTot = - k_*x + omega*t;
  
    // First order contribution
    double Uvert = H/2.0*omega *
                   sinh(k_*(y + h))/sinh(k_*h) *
                   sin(phaseTot);

    // Second order contribution
    Uvert += 3.0/16.0*cel*sq(k_*H)*sinh(2.0*k_*(y + h))
            /pow(sinh(k_*h), 4.0)*sin(2.0*phaseTot);

     return -Uvert;
 }

/*
double wave (double x, double y, double t)
{
  double ak = kWaveData.ak;
  double k_ = kWaveData.k0;
  double h_ = kWaveData.waterdepth;
  double g_ = kWaveData.g;

  double a_ = ak/k_;
  double w = sq(g_*k_*tanh(k_*h_));
  double phaseTot = k_*x - w*t;
  double eta1 = a_*cos(phaseTot);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*phaseTot);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*phaseTot);
  return eta1 + ak*eta2 + sq(ak)*eta3;
}

double u_x (double x, double y, double t)
{
  double ak = kWaveData.ak;
  double k_ = kWaveData.k0;
  double h_ = kWaveData.waterdepth;
  double g_ = kWaveData.g;

  double alpa = 1./tanh(k_*h_);
  double a_ = ak/k_;
  double w = sq(g_*k_*tanh(k_*h_));
  double phaseTot = k_*x - w*t;
  double sgma = sq(g_*k_*tanh(k_*h_)*
		     (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					(sq(alpa) - 1.) + sq(alpa))));
  double A_ = a_*g_/sgma;
  return A_*cosh(k_*(y + h_))/cosh(k_*h_)*k_*cos(phaseTot) +
    ak*3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
    cosh(2.0*k_*(y + h_))*2.*k_*cos(2.0*phaseTot)/cosh(2.0*k_*h_) +
    ak*ak*1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
    (9.*sq(alpa) - 13.)*
    cosh(3.*k_*(y + h_))/cosh(3.*k_*h_)*a_*a_*k_*k_*A_*3.*k_*cos(3.*phaseTot);
}

double u_y (double x, double y, double t)
{
  double ak = kWaveData.ak;
  double k_ = kWaveData.k0;
  double h_ = kWaveData.waterdepth;
  double g_ = kWaveData.g;

  double alpa = 1./tanh(k_*h_);
  double a_ = ak/k_;
  double w = sq(g_*k_*tanh(k_*h_));
  double phaseTot = k_*x - w*t;
  double sgma = sq(g_*k_*tanh(k_*h_)*
		     (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					(sq(alpa) - 1.) + sq(alpa))));
  double A_ = a_*g_/sgma;
  return A_*k_*sinh(k_*(y + h_))/cosh(k_*h_)*sin(phaseTot) +
    ak*3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
    2.*k_*sinh(2.0*k_*(y + h_))*sin(2.0*phaseTot)/cosh(2.0*k_*h_) +
    ak*ak*1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
    (9.*sq(alpa) - 13.)*
    3.*k_*sinh(3.*k_*(y + h_))/cosh(3.*k_*h_)*a_*a_*k_*k_*A_*sin(3.*phaseTot);
}
*/
double u_z (double x, double y, double t)
{
  return 0.0;
}

#endif
