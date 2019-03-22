// Version 20180518-01: Added suggested parameter ranges comments.

// simulation parameters

const double time_step = 1E-4; // time step
const int S = 1E4; // sample interval (1/time_step)
int N_i = 0; // initial timestep
int N_f = 5E6; // final timestep
const int N = N_f - N_i; // number of time steps
const double dx = 1E-2; // lattice spacing in 1st dimension
const int X = 8192; // number of lattice points in 1st dimension
const short int F = 13; // number of fields (phi,h_u,h_d,l,d,8xg)
const int seed = 1; // seed for random numbers
const bool thermal_spectrum = true;

// scales

const double phi0 = 1E5; // 1E8 (phi)
const double l0 = 1E3; // 1E6 (h_u,h_d,l,d)
const double m0 = 1E0; // 1E0 (g)

const double epsilon = 1E-3; // D-term cutoff (m0/l0)

// potential parameters

const double alpha = 0.083; // 0.053 for phi0=1E8 , (1/(log(phi0/m0)+0.5)) in general
const double mphisq = 1.0; // 0.1, 0.01
const double mssq = m0*m0/(phi0*phi0);
const double Tphi_i = 2.0;
const double GammaTphi = 0.01; // 0.1

const double mHusq = -1.0; // 
const double mHdsq = 1.1; // {1.1,...,1.5}
const double mLsq = 0.9; // {0.5,...,0.9}
const double mdsq = 2.0; // {1,...,4}

const double mu = 1.0; // {1,...,4}
const double modAmu = 0.5; // {0,0.5}
const double modAnu = 1.0;
const double modAd = 1.0; // {0,1}

const double argAmu = 0.0; // {0,3.14}
const double argAnu = 0.0;
const double argAd = 0.0; // {0,3.14}

const double ld = 1E2; // ld * (l0/m0) , {1E0,1E1,1E2,1E3}
const double ldg = 0.3;
