#ifndef jet_fns_h
#define jet_fns_h

double circle_area(double th, double dx, double x, double R);
void arange(double array[], int nvals);
void log10spaceWithMidpoints(double *bounds, double *midpoints, double start, double end, int no);
void elecEnergies(double array[], double minEn,double maxEn, int length);
void logspace(double array[], double val_min, double val_max, int n_values);
void gammarange(double array[],double val_min,int n_values);
void gammainterval(double array[], double val_min, double val_max, int n_vals);
void linspace(double array[], double min, double max, double nvals);
void elecErange(double array[],double Emin, double Emax, double nbins);
void cleanpop(double array[], int size_of_array, double cutoff);

int findfirstzero(double array[], int size_of_array);
int findminelement(double array[], int size_of_array);
int findminelementNO0(double array[], int size_of_array);
int findClosest(double *array, double value, int size_array);

double sumno(double value);
double sum_array(double array[], int size_of_array);
double R_s(double M_BH);
float lorentz(float voverc);
double lortovel(double lor);
float theta_lab(float theta_j, float voverc);
float deg2rad(float deg);
float rad2deg(float rad);
float R_new(float R_orig, float theta_open, float x);
float doppler(float voverc, float theta_l);
double freqtoeV(double freq);
double A_PL(double alpha, double E_j, double gamma, double E_min, double E_max, double A_eq);
double R_0(double E_j, double A_eq, double gamma, double B0);
double get_newB(double B0, double R0, double R);
double electron_PL(double A, double alpha, double E_e, double E_max);
double A_new(double A, double Ne_now, double Ne_orig);
double epsilon(double B);
double elec_energy(double eps, double f_c);
double P_total(double B, double beta, double A, double eps, double fq, double alpha, double dx);
double f_crit(double gamma_e, double B);

double emissivity0(double A, double epsilon, double B, double beta, double alpha, double dx, double R);
double emissivity(double A, double epsilon, double B, double beta, double freq,double alpha, double dx, double R);
double emissivity_new(double j0, double freq, double alpha);
double opacity(double j0, double freq, double epsilon, double alpha);
double larmor(double B,double beta, double gamma);
double ICpow(double U, double beta, double gamma);

double new_synchrotron(double R, double dx, double epsilon, double freq, double k);
double new_synchrotron_taylor(double R, double dx, double epsilon, double freq, double k);
double new_synchrotron_ECO(double R, double dx, double epsilon, double freq, double k, double E_e, double E_max);
double new_synchrotron_general(double R, double dx, double epsilon, double freq, double k);
double new_synchrotron_mod(double R, double dx, double epsilon, double freq, double k);
double new_synchrotron_length(double R, double dx, double epsilon, double freq, double k, double x, double L);

double simple_emissivity(double P_tot, double R, double dx);
double j0_from_j(double j, double alpha, double freq);
double k_from_j0(double j0, double freq, double eps, double alpha);
double Iv_from_j_k(double j_v, double k, double R, double dx);
double Iv_from_j_k_small(double j_v, double k, double R, double dx);

double dEe_dfreq(double eps, double freq);
double j_per_hz(double P1, double R, double dNdE, double eps, double freq);
double k_new(double j, double eps, double freq);
double power_emitted(double j, double k, double R);

double maxdx_p(double Ne, double Ee, double eps, double freq, double k, double R);
double maxdx_p_new(double Ne, double Ee, double eps, double freq, double k, double R);
double maxdx_pow(double Ne, double Ee, double eps, double freq, double k, double R, double dx);
double maxdx_p_new2(double Ne, double Ee, double eps, double freq, double k, double R, double E_low);
double dx_P_guess(double E_e, double B, double beta, double gamma);
double dx_P_guess2(double E_e, double B,double beta, double gamma, double E_e2);
double dx_synIC(double P_IC, double P_sync, double beta_bulk, double Ebin_min);
double dx_SIC(double E_e, double Ee2, double B, double beta, double gamma, double U_gam);
double dx_P_simple(double Ne, double E_e, double P_tot);

double fourmom_photon(double E_ph, double phi2);
double dsolid_angle(double angle, double dangle);
double KN_crosssec(double E_ph_p, double phi2);
double KN_longair(double E);
double phden_thin(double L, double freq, double R, double dx);
double phden_thick(double freq, double E_e);

double pluLT(double E, double gamma, double beta,double theta);
double minLT(double E, double gamma, double beta, double theta);
double angLT(double theta, double gamma, double beta);
double invangLT(double theta, double gamma, double beta);
double elecscatt(double E_i, double angle);
double weight_fn(double Ne, double KN, double n_ph, double beta, double theta, double phi2);

double IC_onestep(double beta, double gamma, double Eph,  double theta1, double theta2, double alpha);
double IC_loss(double beta, double gamma, double Urad);

double P_jet(double B, double gamma, double R, double beta);
double P_jet(double B, double gamma, double R, double beta);
double DD_Beffective(double a, double b, double c, double v1, double v2, double v3, double Gamma);
void DD_3Beffective(double a, double b, double c, double v1, double v2, double v3, double Gamma, double *B_effective);
int rand_lim(int limit);
double theta_X(double theta_tot, double theta_circ);
double theta_Y(double theta_tot, double theta_circ);
double Z_perpperp(double cosk, double phik, double Btheta);
double Z_perppara(double cosk, double phik, double Btheta);
double Z_paraperp(double cosk, double phik, double Btheta);
double Z_parapara(double cosk, double phik, double Btheta);
double Z_e(double *_e, double *e, double *k, int perp);
double Sigma_1(double E_elecs,double F_min, double dN_dE, double dEe );
double Sigma_2(double E_elecs,double F_min, double dN_dE, double dEe );


#endif


