#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "jet_fns.h"
#include <assert.h>
#include <gsl/gsl_matrix.h>
#ifdef DEBUG
#define PRINT(message, parameter) \
do { \
    printf(message, parameter); \
} while (0)
#else
#define PRINT(message, parameter)
#endif

#define Me_EV 0.511E6 //electron rest energy
#define C 3.0E8       //speed of light
#define Qe 1.6E-19    //electron charge
#define Me 9.11E-31   //electron mass kg
#define H 6.63E-34    //Plancks constant

/* Program to store some useful Lorentz transform and jets fns */
extern const size_t nLT_sections;
extern const size_t ARRAY_SIZE;


//Calculate effective alpha for a given electron population
void calculate_alpha(double* dN_dE, double* E_elecs, double* effective_alpha){
  for (size_t l = 0; l < ARRAY_SIZE; l++)
  {
      //make alpha dependent on frequency here
      if (l == 0)
      {
          effective_alpha[l] = -(log10(dN_dE[l + 1]) - log10(dN_dE[l])) / (log10(E_elecs[l + 1] * Qe) - log10(E_elecs[l] * Qe));
      }
      else if (l > ARRAY_SIZE - 2)
      {
          effective_alpha[l] = -(log10(dN_dE[l]) - log10(dN_dE[l - 1])) / (log10(E_elecs[l] * Qe) - log10(E_elecs[l - 1] * Qe));
      }
      else
      {
          effective_alpha[l] = -(log10(dN_dE[l + 1]) - log10(dN_dE[l - 1])) / (log10(E_elecs[l + 1] * Qe) - log10(E_elecs[l - 1] * Qe));
      }
      if (isnan(effective_alpha[l]) || isinf(effective_alpha[l]))
      {
          effective_alpha[l] = 8.41283e+01;
      }
  }
}


//Set SSC matrix elements
void set_SSC_matrix(gsl_matrix** A, const size_t N_BLOCKS, double *f_pol,
                    double *f_pol_IC,
                    double *B_effectives,
                    double* E_elecs, double* E_elec_min, double* E_elec_max,
                    double* dEe, double* effective_alpha,
                    double *unitalign0, double *unitalign1, double *unitalign2,
                    double *unitalign3, double *unitalign4, double *unitalign5,
                    double *unitalign6)
{ 

  //cos(th_k) for compton integral
  double cosk[10];
  double phik[10];
  double dcosk = 1 / 5.5;         //1/25.5;
  double dphik = (2 * M_PI) / 10; //(2*M_PI)/ARRAY_SIZE;
  for (int i = 0; i < 10; i++)
  {
      cosk[i] = (i - 5) / 5.5;               //(i-25)/25.5;
      phik[i] = -M_PI + i * (2 * M_PI) / 10; //i*(2*M_PI)/ARRAY_SIZE;
  }
  //parallelize over for loops using openMP
  // #pragma omp parallel for
  for (size_t g = 0; g < N_BLOCKS; g++)
  { //B-fields
      for (size_t h = 0; h < N_BLOCKS; h++)
      {                      
          //now including full RPAR treatment for IC
          int phik_start, phik_end, cosk_start, cosk_end;
          double cosk_single[7];
          double phik_single[7];
          int cosk_list[7];
          int phik_list[7];
          double ang_factor; //reduction in power from small solid angle of zone already accounted for by dx/dxsum in dfactor__, this term remedies
          
          double Btheta = acos(B_effectives[2 + g*3 + h*3*N_BLOCKS + nLT_sections/2*3*N_BLOCKS*N_BLOCKS]);                                       //compton polarization fraction very dependent on this for a single zone, 90deg gives highest
          double zeta = atan(B_effectives[1 + g*3 + h*3*N_BLOCKS + nLT_sections/2*3*N_BLOCKS*N_BLOCKS] / B_effectives[0 + g*3 + h*3*N_BLOCKS + nLT_sections/2*3*N_BLOCKS*N_BLOCKS]); //to rotate each Stokes to lab frame
          // printf("Btheta = %f, zeta = %f\n", Btheta, zeta);

          if (h == g)
          {
              phik_start = 0, phik_end = 10, cosk_start = 0, cosk_end = 10;
              //rotate to B_effective
              ang_factor = 1; //adjust so total IC power always the same
          }
          else if (h != g)
          {
              //first rotate align vector in exactly the same rotation as B -> Beffective, taken care of by unitalign
              //then rotate unitalign vector by -zeta about z axis to get B in x-z plane
              cosk_single[0] = unitalign0[2 + g*3 + h*3*N_BLOCKS];                                                                                                                             //z component isnt affected by -zeta rotation
              phik_single[0] = atan2(unitalign0[0 + g*3 + h*3*N_BLOCKS] * sin(-zeta) + unitalign0[1 + g*3 + h*3*N_BLOCKS] * cos(-zeta), unitalign0[0 + g*3 + h*3*N_BLOCKS] * cos(-zeta) - unitalign0[1 + g*3 + h*3*N_BLOCKS] * sin(-zeta)); //this is between -pi and pi, but phik between 0 and 2pi
              cosk_single[1] = unitalign1[2 + g*3 + h*3*N_BLOCKS];                                                                                                                             //z component isnt affected by -zeta rotation
              phik_single[1] = atan2(unitalign1[0 + g*3 + h*3*N_BLOCKS] * sin(-zeta) + unitalign1[1 + g*3 + h*3*N_BLOCKS] * cos(-zeta), unitalign1[0 + g*3 + h*3*N_BLOCKS] * cos(-zeta) - unitalign1[1 + g*3 + h*3*N_BLOCKS] * sin(-zeta));
              cosk_single[2] = unitalign2[2 + g*3 + h*3*N_BLOCKS]; //z component isnt affected by -zeta rotation
              phik_single[2] = atan2(unitalign2[0 + g*3 + h*3*N_BLOCKS] * sin(-zeta) + unitalign2[1 + g*3 + h*3*N_BLOCKS] * cos(-zeta), unitalign2[0 + g*3 + h*3*N_BLOCKS] * cos(-zeta) - unitalign2[1 + g*3 + h*3*N_BLOCKS] * sin(-zeta));
              cosk_single[3] = unitalign3[2 + g*3 + h*3*N_BLOCKS]; //z component isnt affected by -zeta rotation
              phik_single[3] = atan2(unitalign3[0 + g*3 + h*3*N_BLOCKS] * sin(-zeta) + unitalign3[1 + g*3 + h*3*N_BLOCKS] * cos(-zeta), unitalign3[0 + g*3 + h*3*N_BLOCKS] * cos(-zeta) - unitalign3[1 + g*3 + h*3*N_BLOCKS] * sin(-zeta));
              cosk_single[4] = unitalign4[2 + g*3 + h*3*N_BLOCKS]; //z component isnt affected by -zeta rotation
              phik_single[4] = atan2(unitalign4[0 + g*3 + h*3*N_BLOCKS] * sin(-zeta) + unitalign4[1 + g*3 + h*3*N_BLOCKS] * cos(-zeta), unitalign4[0 + g*3 + h*3*N_BLOCKS] * cos(-zeta) - unitalign4[1 + g*3 + h*3*N_BLOCKS] * sin(-zeta));
              cosk_single[5] = unitalign5[2 + g*3 + h*3*N_BLOCKS]; //z component isnt affected by -zeta rotation
              phik_single[5] = atan2(unitalign5[0 + g*3 + h*3*N_BLOCKS] * sin(-zeta) + unitalign5[1 + g*3 + h*3*N_BLOCKS] * cos(-zeta), unitalign5[0 + g*3 + h*3*N_BLOCKS] * cos(-zeta) - unitalign5[1 + g*3 + h*3*N_BLOCKS] * sin(-zeta));
              cosk_single[6] = unitalign6[2 + g*3 + h*3*N_BLOCKS]; //z component isnt affected by -zeta rotation
              phik_single[6] = atan2(unitalign6[0 + g*3 + h*3*N_BLOCKS] * sin(-zeta) + unitalign6[1 + g*3 + h*3*N_BLOCKS] * cos(-zeta), unitalign6[0 + g*3 + h*3*N_BLOCKS] * cos(-zeta) - unitalign6[1 + g*3 + h*3*N_BLOCKS] * sin(-zeta));

              for (size_t n = 0; n < 7; n++)
              {
                  cosk_list[n] = findClosest(cosk, cosk_single[n], 10);
                  phik_list[n] = findClosest(phik, phik_single[n], 10);
              }
              cosk_start = 0, cosk_end = 7;
              phik_start = 0, phik_end = 1;
              ang_factor = 14.3; //ie 20*5 = 100 = size(phik) * size(cosk)
          }
          else
          {
              continue;
          }
          
          #pragma omp parallel for collapse(4)
          for (int n = 0; n < ARRAY_SIZE; n++)
          { //f_polIC (compton energy)
              for (int l = 0; l < ARRAY_SIZE; l++)
              { //f_pol (sync energy)
                  for (int p = phik_start; p < phik_end; p++)
                  { //phi
                      for (int m = cosk_start; m < cosk_end; m++)
                      { //cosk
                          int nn = m;
                          int o;
                          if (h != g)
                          {
                              p = phik_list[nn];
                              m = cosk_list[nn];
                          }

                          double F_min = sqrt(f_pol_IC[n] / (2 * f_pol[l] * (1 - cosk[m])));
                          double v_k[3];     //incoming photon vector
                          double e_k[3];     //incoming polarization vector perpendicular to B
                          double _e[2];      //outgoing polarization vector
                          double e_kpara[3]; //incoming polarization vector parallel to B

                          v_k[0] = sqrt(1 - cosk[m] * cosk[m]) * cos(phik[p]), v_k[1] = sqrt(1 - cosk[m] * cosk[m]) * sin(phik[p]), v_k[2] = cosk[m]; //incoming photon direction vector
                          e_k[0] = sqrt(1 - cosk[m] * cosk[m]) * sin(phik[p]) * cos(Btheta);
                          e_k[1] = sin(Btheta) * cosk[m] - cos(Btheta) * sqrt(1 - cosk[m] * cosk[m]) * cos(phik[p]); //incoming photon polarization vector (perp)
                          e_k[2] = -sin(Btheta) * sqrt(1 - cosk[m] * cosk[m]) * sin(phik[p]);
                          //this is just cross product of e_k and v_k above
                          e_kpara[0] = pow(cosk[m], 2) * sin(Btheta) - cos(Btheta) * cos(phik[p]) * cosk[m] * sqrt(1 - cosk[m] * cosk[m]) + (1 - cosk[m] * cosk[m]) * pow(sin(phik[p]), 2) * sin(Btheta),
                          e_kpara[1] = -sin(Btheta) * (1 - cosk[m] * cosk[m]) * sin(2 * phik[p]) / 2 - cosk[m] * sqrt(1 - cosk[m] * cosk[m]) * sin(phik[p]) * cos(Btheta), //incoming photon polarization vector (para)
                          e_kpara[2] = (1 - cosk[m] * cosk[m]) * pow(sin(phik[p]), 2) * cos(Btheta) - sqrt(1 - cosk[m] * cosk[m]) * cos(phik[p]) * sin(Btheta) * cosk[m] + (1 - cosk[m] * cosk[m]) * pow(cos(phik[p]), 2) * cos(Btheta);
                          //e_k above is not normalised
                          double q_theta = pow(1 - pow(cos(Btheta) * cosk[m] + sin(Btheta) * cos(phik[p]) * sqrt(1 - cosk[m] * cosk[m]), 2), (effective_alpha[l] + 1) / 4);

                          if (F_min * (Me_EV) > E_elecs[ARRAY_SIZE - 1])
                          {}

                          else if (F_min * (Me_EV) > E_elecs[0])
                          {
                              for (o = 0; o < ARRAY_SIZE; o++)
                              { //find Fmin location in E_elecs
                                  if (F_min * (Me_EV) > E_elec_min[o] && F_min * (Me_EV) < E_elec_max[o])
                                  {
                                      break;
                                  }
                              }

                              for (int i = o; i < ARRAY_SIZE; i++)
                              { //E_e
                                  //initial constant factor to make sure this matches with isotropic electron losses, original constant from paper is 1.2859E-91

                                  double Pperpperp = 3.95417e-103 * q_theta * /*dfreqs_polIC[n] */ f_pol_IC[n] / f_pol[l] * F_min * /*dfreqs_pol[l]* */ (1 / (f_pol[l] * H)) * dcosk * dphik * ang_factor //this is Power(freq), multiply by freq later to get nuF(nu)
                                              * (Z_e(_e, e_k, v_k, 1) * (Sigma_2(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i]) + Sigma_1(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i])) + Sigma_2(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i]));
                                  //printf("Pperperp %.5e\n", Pperpperp);
                                  gsl_matrix_set(A[1+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1, 
                                                  gsl_matrix_get(A[1+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1) + (Pperpperp / (N_BLOCKS)) * cos(2 * atan2(_e[1] * cos(zeta) + _e[0] * sin(zeta), _e[0] * cos(zeta) - _e[1] * sin(zeta))));
                                  gsl_matrix_set(A[2+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1, 
                                                  gsl_matrix_get(A[2+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1) + (Pperpperp / (N_BLOCKS)) * sin(2 * atan2(_e[1] * cos(zeta) + _e[0] * sin(zeta), _e[0] * cos(zeta) - _e[1] * sin(zeta))));
                                                                    
                                  double Pperppara = 3.95417e-103 /*8.26784E-118*/ * q_theta * /*dfreqs_polIC[n] */ f_pol_IC[n] / f_pol[l] * F_min * /*dfreqs_pol[l]* */ (1 / (f_pol[l] * H)) * dcosk * dphik * ang_factor * (Z_e(_e, e_kpara, v_k, 1) * (Sigma_2(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i]) + Sigma_1(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i])) + Sigma_2(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i]));

                                  gsl_matrix_set(A[1+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2, 
                                                  gsl_matrix_get(A[1+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2)+ (Pperppara / (N_BLOCKS)) * cos(2 * atan2(_e[1] * cos(zeta) + _e[0] * sin(zeta), _e[0] * cos(zeta) - _e[1] * sin(zeta))));
                                  gsl_matrix_set(A[2+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2, 
                                                  gsl_matrix_get(A[2+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2) + (Pperppara / (N_BLOCKS)) * sin(2 * atan2(_e[1] * cos(zeta) + _e[0] * sin(zeta), _e[0] * cos(zeta) - _e[1] * sin(zeta))));

                                  double Pparaperp = 3.95417e-103 /*8.26784E-118*/ * q_theta * /*dfreqs_polIC[n] */ f_pol_IC[n] / f_pol[l] * F_min * /*dfreqs_pol[l]* */ (1 / (f_pol[l] * H)) * dcosk * dphik * ang_factor * (Z_e(_e, e_k, v_k, 0) * (Sigma_2(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i]) + Sigma_1(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i])) + Sigma_2(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i]));
                                  gsl_matrix_set(A[1+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1, 
                                                  gsl_matrix_get(A[1+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1) + (Pparaperp / (N_BLOCKS)) * cos(2 * atan2(_e[1] * cos(zeta) + _e[0] * sin(zeta), _e[0] * cos(zeta) - _e[1] * sin(zeta)))); 
                                  gsl_matrix_set(A[2+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1, 
                                                  gsl_matrix_get(A[2+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1) + (Pparaperp / (N_BLOCKS)) * sin(2 * atan2(_e[1] * cos(zeta) + _e[0] * sin(zeta), _e[0] * cos(zeta) - _e[1] * sin(zeta)))); 
                                                  
                                  double Pparapara = 3.95417e-103 /*8.26784E-118*/ * q_theta * /*dfreqs_polIC[n] */ f_pol_IC[n] / f_pol[l] * F_min * /*dfreqs_pol[l]* */ (1 / (f_pol[l] * H)) * dcosk * dphik * ang_factor * (Z_e(_e, e_kpara, v_k, 0) * (Sigma_2(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i]) + Sigma_1(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i])) + Sigma_2(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i]));
                                  gsl_matrix_set(A[1+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2, 
                                                  gsl_matrix_get(A[1+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2) + (Pparapara / (N_BLOCKS)) * cos(2 * atan2(_e[1] * cos(zeta) + _e[0] * sin(zeta), _e[0] * cos(zeta) - _e[1] * sin(zeta)))); 
                                  gsl_matrix_set(A[2+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2, 
                                                  gsl_matrix_get(A[2+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2) + (Pparapara / (N_BLOCKS)) * sin(2 * atan2(_e[1] * cos(zeta) + _e[0] * sin(zeta), _e[0] * cos(zeta) - _e[1] * sin(zeta)))); 

                                  gsl_matrix_set(A[0+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2, 
                                                  gsl_matrix_get(A[0+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2) + (Pparapara / (N_BLOCKS)) + (Pperppara / (N_BLOCKS)));
                                  gsl_matrix_set(A[0+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1, 
                                                  gsl_matrix_get(A[0+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1) + (Pparaperp / (N_BLOCKS)) + (Pperpperp / (N_BLOCKS)));
 
                              }
                          }

                          else
                          {   

                              for (int i = 0; i < ARRAY_SIZE; i++)
                              { //E_e

                                  double Pperpperp = 3.95417e-103 * q_theta * /*dfreqs_polIC[n] */ f_pol_IC[n] / f_pol[l] * F_min * /*dfreqs_pol[l]* */ (1 / (f_pol[l] * H)) * dcosk * dphik * ang_factor //this is Power(freq), multiply by freq later to get nuF(nu)
                                              * (Z_e(_e, e_k, v_k, 1) * (Sigma_2(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i]) + Sigma_1(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i])) + Sigma_2(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i]));
                                  //printf("Pperperp %.5e\n", Pperpperp);
                                  gsl_matrix_set(A[1+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1, 
                                                  gsl_matrix_get(A[1+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1) + (Pperpperp / (N_BLOCKS)) * cos(2 * atan2(_e[1] * cos(zeta) + _e[0] * sin(zeta), _e[0] * cos(zeta) - _e[1] * sin(zeta))));
                                  gsl_matrix_set(A[2+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1, 
                                                  gsl_matrix_get(A[2+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1) + (Pperpperp / (N_BLOCKS)) * sin(2 * atan2(_e[1] * cos(zeta) + _e[0] * sin(zeta), _e[0] * cos(zeta) - _e[1] * sin(zeta))));
                                                                    
                                  double Pperppara = 3.95417e-103 /*8.26784E-118*/ * q_theta * /*dfreqs_polIC[n] */ f_pol_IC[n] / f_pol[l] * F_min * /*dfreqs_pol[l]* */ (1 / (f_pol[l] * H)) * dcosk * dphik * ang_factor * (Z_e(_e, e_kpara, v_k, 1) * (Sigma_2(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i]) + Sigma_1(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i])) + Sigma_2(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i]));

                                  gsl_matrix_set(A[1+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2, 
                                                  gsl_matrix_get(A[1+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2)+ (Pperppara / (N_BLOCKS)) * cos(2 * atan2(_e[1] * cos(zeta) + _e[0] * sin(zeta), _e[0] * cos(zeta) - _e[1] * sin(zeta))));
                                  gsl_matrix_set(A[2+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2, 
                                                  gsl_matrix_get(A[2+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2) + (Pperppara / (N_BLOCKS)) * sin(2 * atan2(_e[1] * cos(zeta) + _e[0] * sin(zeta), _e[0] * cos(zeta) - _e[1] * sin(zeta))));

                                  double Pparaperp = 3.95417e-103 /*8.26784E-118*/ * q_theta * /*dfreqs_polIC[n] */ f_pol_IC[n] / f_pol[l] * F_min * /*dfreqs_pol[l]* */ (1 / (f_pol[l] * H)) * dcosk * dphik * ang_factor * (Z_e(_e, e_k, v_k, 0) * (Sigma_2(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i]) + Sigma_1(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i])) + Sigma_2(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i]));
                                  gsl_matrix_set(A[1+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1, 
                                                  gsl_matrix_get(A[1+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1) + (Pparaperp / (N_BLOCKS)) * cos(2 * atan2(_e[1] * cos(zeta) + _e[0] * sin(zeta), _e[0] * cos(zeta) - _e[1] * sin(zeta)))); 
                                  gsl_matrix_set(A[2+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1, 
                                                  gsl_matrix_get(A[2+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1) + (Pparaperp / (N_BLOCKS)) * sin(2 * atan2(_e[1] * cos(zeta) + _e[0] * sin(zeta), _e[0] * cos(zeta) - _e[1] * sin(zeta)))); 
                                                  
                                  double Pparapara = 3.95417e-103 /*8.26784E-118*/ * q_theta * /*dfreqs_polIC[n] */ f_pol_IC[n] / f_pol[l] * F_min * /*dfreqs_pol[l]* */ (1 / (f_pol[l] * H)) * dcosk * dphik * ang_factor * (Z_e(_e, e_kpara, v_k, 0) * (Sigma_2(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i]) + Sigma_1(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i])) + Sigma_2(E_elecs[i], Me_EV * Qe * F_min, 1, dEe[i]));
                                  gsl_matrix_set(A[1+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2, 
                                                  gsl_matrix_get(A[1+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2) + (Pparapara / (N_BLOCKS)) * cos(2 * atan2(_e[1] * cos(zeta) + _e[0] * sin(zeta), _e[0] * cos(zeta) - _e[1] * sin(zeta)))); 
                                  gsl_matrix_set(A[2+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2, 
                                                  gsl_matrix_get(A[2+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2) + (Pparapara / (N_BLOCKS)) * sin(2 * atan2(_e[1] * cos(zeta) + _e[0] * sin(zeta), _e[0] * cos(zeta) - _e[1] * sin(zeta)))); 

                                  gsl_matrix_set(A[0+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2, 
                                                  gsl_matrix_get(A[0+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2) + (Pparapara / (N_BLOCKS)) + (Pperppara / (N_BLOCKS)));
                                  gsl_matrix_set(A[0+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1, 
                                                  gsl_matrix_get(A[0+h*3], n, (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1) + (Pparaperp / (N_BLOCKS)) + (Pperpperp / (N_BLOCKS)));
                              }
                          }

                          if (h != g)
                          {
                              m = nn;
                          }
                      }
                  }
              }
          }
      }
  }
}



//Calculates the area of intersection between circle of radius R (jet radius) and circle of radius x
double circle_area(double d, double dx, double x, double R){
    // d is the distance between centre points of the two circles
    if (x <= R) {
       double area;
       double r = x;
       double r_low = x-dx; //required for annulus calculation

       if (d <= R - r){ //whole emission circle fits in jet circle
           if (x == dx) { //first contribution is full circle
               return area = M_PI * pow(dx,2);
           } else {
               return area = M_PI * (pow(x,2) - pow((x-dx),2)); //area of annulus
           }
       }
       else { //emission circle overlaps jet circle, should always be annulus in this case
           double d_1 = (pow(R,2) - pow(r,2) + pow(d,2)) / (2 * d);
           double d_2 = d - d_1;
           double area_low;
           double d_2low;
           double d_1low;

           if (d <= R - r_low) {//case where inner circle is within jet but outer peeks out
             area_low = M_PI * pow(r_low,2);
           } else {
             d_1low = (pow(R,2) - pow(r_low,2) + pow(d,2)) / (2 * d);
             d_2low = d - d_1low;

             area_low = pow(R,2) * acos(d_1low/R) - d_1low * sqrt(pow(R,2) - pow(d_1low,2)) + pow(r_low,2) * acos(d_2low/r_low) - d_2low * sqrt(pow(r_low,2) - pow(d_2low,2));
           }

           return area = pow(R,2) * acos(d_1/R) - d_1 * sqrt(pow(R,2) - pow(d_1,2)) + pow(r,2) * acos(d_2/r) - d_2 * sqrt(pow(r,2) - pow(d_2,2)) - area_low; //annulus area
       }
    } else {
	      double area;
        double r = R;
        double Rbig = x; //required for annulus calculation
	      double Rbig_low = x-dx;

	      if (d <= Rbig - r){ //whole emission circle fits in jet circle
	        printf("radius of emission bigger than jet radius: ERROR");
          return 0;
        }
        else { //emission circle overlaps jet circle, should always be annulus in this case
           double d_1 = (pow(Rbig,2) - pow(r,2) + pow(d,2)) / (2 * d);
           double d_2 = d - d_1;
           double area_low;
           double d_2low;
           double d_1low;


           d_1low = (pow(Rbig_low,2) - pow(r,2) + pow(d,2)) / (2 * d);
           d_2low = d - d_1low;

           area_low = pow(Rbig_low,2) * acos(d_1low/Rbig_low) - d_1low * sqrt(pow(Rbig_low,2) - pow(d_1low,2)) + pow(r,2) * acos(d_2low/r) - d_2low * sqrt(pow(r,2) - pow(d_2low,2));
          

           return area = pow(Rbig,2) * acos(d_1/Rbig) - d_1 * sqrt(pow(Rbig,2) - pow(d_1,2)) + pow(r,2) * acos(d_2/r) - d_2 * sqrt(pow(r,2) - pow(d_2,2)) - area_low; //annulus area
        }
	
    }

}

void arange(double array[], int nvals)
/* Analogy to numpy.linspace */
{
  double element = 0.0;
  int i;
  for (i=0; i<nvals; i++)
    {
    array[i] = element;
    element += 1.0;
    }
}


/* Fn to return a logarithmically spaced array of bounds and midpoint
* bounds size is one greater than midpoints!
* start < end
*/
void log10spaceWithMidpoints(double *bounds, double *midpoints, double start, double end, int no)
{
  int j;

  bounds[0] = start;
  //printf("no %.2d", no);
  for (j=1; j<no; j++)//element<=end; j++)
    {
      bounds[j] = pow(10,log10(start) + j*(log10(end)-log10(start))/(no-1));
      midpoints[j-1] = pow(10, log10(bounds[j-1]) + (log10(bounds[j])-log10(bounds[j-1]))/2);

   }

}


/* Fn to define an array range of energies */
void elecEnergies(double array[], double minEn, double maxEn, int length)
{
  double stepsize, element;
  int i, arraysize;
  stepsize = (maxEn-minEn)/(length-1);
  //  arraysize = sizeof(array);//sizeof(*array);//sizeof(array[0]);
  //  printf("\nThe stepsize is: %.2e", stepsize);
  //printf("\n array size is: %d", arraysize);
  element = minEn; //initialize
    for (i=0; i<=length; i++)
      {
	array[i] = element;
	element += stepsize;
      }
}


/*Fn to return linearly spaced numbers */
void linspace(double array[], double minval, double maxval, double nvals)
{
  double element_to_add, stepsize;
  int j; //loop parameter
  stepsize = (maxval-minval)/nvals;
  element_to_add = minval;

  for (j=0; j<nvals; j++)
    {
      array[j] = element_to_add;
      element_to_add += stepsize;
    }
}


void elecErange(double array[], double Emin, double Emax, double nbins)
{
  double element, ratio;
  int i;
  //  ratio = (1.3/0.7); // defined like this in paper
  ratio = pow(Emax/Emin, 1/nbins);
  //printf("ratio is: %.5e \n", ratio);
  element = Emin;
  for (i=0; i<nbins; i++)
    {
      //printf("%.5e \n", element);
      array[i] = element;
      //printf("%.5e \n", array[i]);
      element*=ratio;
    }
}


/* Fn to logarithmically space n values between val_min & val_max */
void logspace(double array[], double val_min, double val_max, int n_values)
{
  double ratio, element, n_val_float;
  int i;
  n_val_float = n_values;
  ratio = pow((val_max/val_min), (1/n_val_float)); //gets the ratio to multiply by
  element = val_min; //initialize
  //printf("\n ratio, element= %.2e\t%.2e \n", ratio, element);
  //printf("\n max/min: %.2e\t%.2e\t%.2e\t%.2d\n", val_max, val_min, val_max/val_min, 1/n_values);///val_min);

  for (i=0; i<=n_values; i++)//element<val_max*ratio; i++)//i<n_values; i++)
    {
      array[i] = element;
      element *= ratio;
    }
}


/* Fn to generate spaced gamma values */
void gammarange(double array[], double val_min, int n_values)
{
  int i;
  double element;
  element = val_min;
  for (i=0; i<n_values; i++)
    {
      array[i] = element;
      element*=(1.3/0.7);
    }
}

/* As above, but specify max and min values */
void gammainterval(double array[], double val_min, double val_max, int n_vals)
{
  int i;
  double element, interval;
  interval = (val_max-val_min)/n_vals;
  element = val_min;
  for (i=0; i<n_vals; i++)
    {
      array[i] = element;
      element += interval;
    }
}


/* Fn to remove any electron populations to avoid nans in code */
void cleanpop(double array[], int size_of_array, double cutoff) //cutoff sets the upper limit to going to zero
{
  int i; //loop parameter
  for (i=0; i<size_of_array; i++)
    {
      if (array[i]<cutoff)
	{
	  array[i] = 0.0; //clean up
	}
    }
}


/* Fn to find the first element with value 0, or the lowest element */
int findfirstzero(double array[], int size_of_array)
{
  int i, element=size_of_array-1; //find the first zero element
  double min=1E300; //initialise as some ridiculously huge value
  for (i=size_of_array-1; i>=0; i--)//no -1 before
    {
      if (array[i]==0)
	{
	  element = i;
	  //printf("%d \n", element);
	}
    }
  if (element==size_of_array-1)//find the minimum value instead +1 before
    {
      for (i=0; i<size_of_array; i++) //no -1 before
	{
	  if (array[i]<min)
	    {
              //printf("%.2e\t%s\t%.2e\t%d \n", array[i], "<", min, element);	     
              min = array[i];
	      element = i;
	    }
	}
    }
  return element;
}


/*Fn to find the minimum element in an array */
int findminelement(double array[], int size_of_array)
{
  int i, element=size_of_array-1; //need -1 as array size x labelled 0->x-1
  double min = 1E300; //initialise BIG

  for (i=0; i<size_of_array; i++)
    {
      if (array[i]<min) ///&& array[i]!=0.0) //needs to be non zero or code won't work!
	{
	  //printf("%.2e\t%s\t%.2e\t%d \n", array[i], "<", min, element);
	  min = array[i];
	  element = i;
	}
    }
  return element;
}

int findminelementNO0(double array[], int size_of_array)
{
  int i, element=size_of_array-1; //need -1 as array size x labelled 0->x-1
  double min = 1E300; //initialise BIG

  for (i=0; i<size_of_array; i++){
        if (array[i]<min && array[i]!=0.0) //needs to be non zero or code won't work!
        {
        //printf("%.2e\t%s\t%.2e\t%d \n", array[i], "<", min, element);
        min = array[i];
        element = i;
        }
    }
  return element;
}

int findClosest(double *array, double value, int size_array){ //finds element of member of array which is closest to value
  int i;
  double func_array[size_array];
  for (i=0; i<size_array; i++){
    func_array[i] = fabs(array[i] - value);
  }
  return findminelement(func_array, size_array);
}


/*Fn to compute the sum of numbers */
double sumno(double value)
{
  double sum=0;
  int i;

  for (i=0; i<=value; i++)
    {
      sum+=i;
    }
  return sum;
}


/*Fn to find the sum of all elemets in an array */
double sum_array(double array[], int size_of_array)
{
  int i;
  double sum=0;

  for (i=0; i<size_of_array; i++)
    {
      sum+=array[i];
    }
  return sum;
}


/* Fn to convert degrees to radians */
float deg2rad(float deg)
{
  float rad;

  rad = deg*M_PI/180.0;
  return (rad);
}

/*Fn to convert radians to degress */
float rad2deg(float rad)
{
  float deg;
  deg = rad*180/M_PI;
  return (deg);
}

/* Fn for the Lorentz factor */
float lorentz(float voverc)
{
  float ans;
  ans = pow((1.0-pow(voverc, 2.0)), -0.5);
  return (ans);
}


/*Fn to return velocity from Lorentz factor */
double lortovel(double lor)
{
  double ans;
  assert (lor>=1);
  ans = pow((1.0-1.0/(lor*lor)), 0.5);
  return (ans);
}


/*Fn to find the Scwarzchild radius */
double R_s(double M_BH)
{
  double ans;
  ans = (2.0*6.67E-11*M_BH)/(3.0E8*3.0E8);
  return ans;
}


/*Fn to propagate the jet radius */
float R_new(float R_orig, float theta_open, float x)
{
  float ans;
  ans = R_orig + x*tan(theta_open);
  return ans;
}

/* Fn for the angle transform */
float theta_lab(float theta_j, float voverc)
{
  float ans;
  ans = asin(sin(theta_j)/(lorentz(voverc)*(1+voverc*cos(theta_j))));

  /* apply correction for oppositely propagating jet */
  if (fabsf(theta_j) > 90.0*M_PI/180.0)
    {
    ans = M_PI-ans;
    }

  return (ans);
}


/* Fn to define the Doppler factor (INPUT IN RADIANS) */
float doppler(float voverc, float theta_l)//, char ang="r")
{
  float ans;
  ans = 1.0/(lorentz(voverc)*(1.0-voverc*cos(theta_l)));
  return (ans);
}


/*Fn to convert frequency to eV */ 
double freqtoeV(double freq)
{
  double ans;
  ans = 6.63E-34*freq/1.6E-19;
  return (ans);
}


/* THE BELOW SECTION HAS THE AIM OF REPRODUCING SYNCHROTRON SPECTRA */ 

/*Fn to return the A PL coefficient GAMMA IS GAMMA BULK */
double A_PL(double alpha, double E_j, double gamma, double E_min, double E_max, double A_eq)
{
  double ans;
  ans = ((2.0-alpha)*E_j)/(gamma*gamma*(1.0+A_eq)*(pow(E_max, (2.0-alpha))-pow(E_min, (2.0-alpha))));
  return (ans);
}


/* Function to determine jet radius at the base */
double R_0(double E_j, double A_eq, double gamma, double B0)
{
  /* need to define mu0 and lc (light second) in main code */
  double ans, sqrt;
  sqrt = (2.0*E_j*A_eq*4.0*M_PI*pow(10.0, -7.0))/(gamma*gamma*(M_PI*3.0E8*B0*B0)*(1.0+A_eq));
  ans = pow(sqrt, 0.5);
  return (ans);
}


/* Fn to scale the magnetic field as the jet radius increases */
double get_newB(double B0, double R0, double R)
{
  double ans;
  ans = B0*(R0/R);
  return ans;
}


/* Fn to define the PL energy distribution of the electron population */
double electron_PL(double A, double alpha, double E_e, double E_max)
{
  double ans;
  ans = A*pow(E_e, -alpha)*pow(M_E, -(E_e/E_max));
  return ans;
}


/* Fn to scale A as the jet radius increases */
double A_new(double A, double Ne_now, double Ne_orig)
{
  double ans;
  ans = A*(Ne_now/Ne_orig); /* Ne_orig defined at x=0 */
  return (ans);
}


/* Fn to define the critical frequency at which all electrons emit at */
double f_crit(double gamma_e, double B)
{
  double ans;
  ans = (3.0*gamma_e*gamma_e*1.6E-19*B)/(4.0*M_PI*9.11E-31);
  return (ans);
}


/* Fn to define epsilon (eqn. 2.18b in PC paper a) */
double epsilon(double B)
{
  double ans;
  ans = (4.0*M_PI*pow(9.11E-31, 3.0)*pow(3.0E8, 4.0))/(3.0*1.6E-19*B);
  return ans;
}


/* Fn to calculate the energy of the electron population */
double elec_energy(double eps, double f_c)
{
  double ans;
  ans = pow((eps*f_c), 0.5);
  return (ans);
}


/* Fn to define total power in section dx as a fn of frequency */
double P_total(double B, double beta, double A, double eps, double fq, double alpha, double dx)
{
  double ans;
  ans = (6.65E-29*B*B*beta*beta*A*eps*pow((eps*fq), (1.0-alpha)/2.0)*dx)/(3*4*M_PI*10E-7*pow(9.11E-31, 2)*pow(3E8, 3)*3E8);
  return (ans);
}



/* Now include some functions to calculate the opacity through jet sections */

/* Fn to get the emissivity coefficient */
double emissivity0(double A, double epsilon, double B, double beta, double alpha, double dx, double R)
{
  double ans;
  ans = (A*epsilon*6.65E-29*B*B*beta*beta*pow(epsilon, (1-alpha)/2))/(3*M_PI*R*R*4*M_PI*1E-7*pow(9.11E-31, 2)*pow(3E8, 4.0));
  return ans;
}


/* Fn to get j0 times freq */
double emissivity(double A, double epsilon, double B, double beta, double freq, double alpha, double dx, double R)
{
  double ans;
  ans = (A*epsilon*6.65E-29*B*B*beta*beta*pow((epsilon*freq), ((1-alpha)/2)))/(3*M_PI*R*R*4*M_PI*1E-7*pow(9.11E-31, 2.0)*pow(3E8, 4.0)); //c^4 as lc=c metres.
  return ans;
}

/* Simplified above */
double emissivity_new(double j0, double freq, double alpha)
{
  double ans;
  ans = j0*pow(freq, (1.0-alpha)/2.0);
  return ans;
}
     						

/* Fn to calculate the frequency dependent opacity of a jet section */
double opacity(double j0, double freq, double epsilon, double alpha)
{
  double ans;
  ans = (j0*3E8*3E8*pow(freq, (-alpha-4.0)/2.0))/(2*pow(epsilon, 0.5));
  return ans;
}	 


/*Fn to calculate the sync power from one electron */
double larmor(double B, double beta, double gamma)
{
  double ans;
  ans = 2.0*6.65E-29*3.0E8*B*B*beta*beta*gamma*gamma/(3.0*4.0*M_PI*1.0E-7);
  return ans;
}


/*Fn to calculate IC power */
double ICpow(double U, double beta, double gamma)
{
  double ans;
  ans = (4.0/3.0)*6.65E-29*3.0E8*beta*beta*gamma*gamma*U;
  return ans;
}


/*Fn to calculate the original synchrotron */
double old_syn(double P_single, double N_e, double dx)
{
  double ans;
  ans = P_single*N_e*dx/3E8;
  return ans;
}


/*Fn to calculate the opacity modified synchrotron */
double new_synchrotron(double R, double dx, double epsilon, double freq, double k)
{
  double ans;
  ans = M_PI*R*dx*(2*pow(epsilon, 0.5))*pow(freq, (5.0/2.0))*(1.000-pow(M_E, -k*R))/(3E8*3E8);
  return ans;
}


/*Fn to calculate the opacity modified synchrotron */
double new_synchrotron_taylor(double R, double dx, double epsilon, double freq, double k)
{
  double ans;
  ans = M_PI*R*dx*((2*pow(epsilon, 0.5)*pow(freq, (5.0/2.0)))/(3E8*3E8))*(k*R);
  return ans;
}



/*Fn to calculate the opacity modified synchrotron with an ECO */
double new_synchrotron_ECO(double R, double dx, double epsilon, double freq, double k, double E_e, double E_max)
{
  double ans, param;
  param = (1.000-pow(M_E, -k*R));

  if (param == 0.0)
    {
      ans =  pow(M_E, (-E_e/E_max))*M_PI*R*dx*(2*pow(epsilon, 0.5)*pow(freq, (5.0/2.0)))*(k*R)/(3E8*3E8);
    }
  else
    {
      ans = pow(M_E, (-E_e/E_max))*M_PI*R*dx*((2*pow(epsilon, 0.5)*pow(freq, (5.0/2.0)))/(3E8*3E8))*(1.000-pow(M_E, -k*R));
    }
  return ans;
}


/*Fn to get the synchrotron at any time */
/*Fn to calculate the opacity modified synchrotron */
double new_synchrotron_general(double R, double dx, double epsilon, double freq, double k)
{
  double ans, param;
  param = (1.000-pow(M_E, -k*R));

  if (param == 0.0)
    {
      ans = M_PI*R*dx*(2*pow(epsilon, 0.5)*pow(freq, (5.0/2.0)))*(k*R)/(3E8*3E8);
    }
  else
    {
      ans = M_PI*R*dx*2*pow(epsilon, 0.5)*pow(freq, (5.0/2.0))*(1.000-pow(M_E, -k*R))/(3E8*3E8);
    } 
 return ans;
}


/* eqn 2.28 in thesis */
double new_synchrotron_mod(double R, double dx, double epsilon, double freq, double k)
{
  double ans, param;
  param = (1.000-pow(M_E, -k*R));

  if (param == 0.0)
    {
      ans = M_PI*R*R*(2*pow(epsilon, 0.5)*pow(freq, (5.0/2.0)))*(k*dx)/(3E8*3E8);
    }
  else
    {
      ans = M_PI*R*R*2*pow(epsilon, 0.5)*pow(freq, (5.0/2.0))*(1.000-pow(M_E, -k*dx))/(3E8*3E8);
    }
  return ans;
}

/*Fn to get the synchrotron at any time */
/*Fn to calculate the opacity modified synchrotron */
double new_synchrotron_length(double R, double dx, double epsilon, double freq, double k, double x, double L)
{
  double ans, param;
  param = (1.000-pow(M_E, -k*(L-x)));

  if (param == 0.0)
    {
      ans = M_PI*R*dx*(2*pow(epsilon, 0.5)*pow(freq, (5.0/2.0)))*(k*(L-x))/(3E8*3E8);
    }
  else
    {
      ans = M_PI*R*dx*2*pow(epsilon, 0.5)*pow(freq, (5.0/2.0))*(1.000-pow(M_E, -k*(L-x)))/(3E8*3E8);
    }
  return ans;
}




/*Some simplified functions for the new code*/

/*fn to calculate total emissivity */
double simple_emissivity(double P_tot, double R, double dx)
{
  double ans;
  ans = P_tot/(M_PI*R*R*dx);
  return ans;
}


/* Fn to get j0 from j */
double j0_from_j(double j, double alpha, double freq)
{
  double ans;
  ans = j/(pow(freq, ((1.0-alpha)/2.0)));
  return ans;
}

/*Fn to get opacity from j0*/
double k_from_j0(double j0, double freq, double eps, double alpha)
{
  double ans;
  ans = j0*3E8*3E8*pow(freq, -1*(alpha+4.0)/2.0)/(2*pow(eps, 0.5));
  return ans;
}

/* Fn to get intensity using emissivity and opacity */
double Iv_from_j_k(double j_v, double k, double R, double dx)
{
  double ans;
  ans = M_PI*R*dx*(j_v/k)*(1-pow(M_E, -k*R));
  return ans;
} 

/* Fn to get intensity using emissivity and opacity */
double Iv_from_j_k_small(double j_v, double k, double R, double dx)
{
  double ans;
  ans = M_PI*R*dx*(j_v/k)*(k*R);
  return ans;
}



/*Some new fns to redo synchrotron properly */

/* Fn to get the electron energy per frequency bin */
double dEe_dfreq(double eps, double freq)
{
  double ans;
  ans = 0.5*pow((eps/freq), 0.50);
  return ans;
}


/* Fn to get the emissivity per unit volume per hertz */
double j_per_hz(double P1, double R, double dNdE, double eps, double freq)
{
  double ans;
  ans = P1*(1.0/(M_PI*R*R*3.0E8))*dNdE*0.5*pow((eps/freq), 0.50);
  return ans;
}


/* Fn to defone opacity based on emissivity */
double k_new(double j, double eps, double freq)
{
  double ans;
  ans = (j*3E8*3E8)/(2.0*pow(eps, 0.50)*pow(freq, 2.50));
  return ans;
}


/* Fn to compute observed emission */
double power_emitted(double j, double k, double R)
{
  double ans, param;
  param = (1.000-pow(M_E, -k*R));
  if (param == 0.0)
    {
      ans = (j*R);
    }
  else
    {
      ans = (j/k)*(1-pow(M_E, -k*R));
      //printf("Using (b) \n");
    }
  return ans;
}


/*Fn to get max allowed dx from power constraints */
double maxdx_p(double Ne, double Ee, double eps, double freq, double k, double R)
{
  double ans, sqrt, param;
  param = (1.000-pow(M_E, -k*R));

  if (param == 0.0)
    {
      sqrt = (Ne*pow(3E8, 3.0)*Ee)/(2.0*M_PI*R*pow(eps, 0.5)*pow(freq, 2.5)*(k*R));
    }
  else
    {
      sqrt = (Ne*pow(3E8, 3.0)*Ee)/(2.0*M_PI*R*pow(eps, 0.5)*pow(freq, 2.5)*(1.00-pow(M_E, -k*R)));
    }
  ans = pow(sqrt, 0.5);
  return ans;
}


/*Fn to get max allowed dx from power constraints */
double maxdx_p_new(double Ne, double Ee, double eps, double freq, double k, double R)
{
  double sqrt, param;
  param = (1.000-pow(M_E, -k*R));

  if (param == 0.0)
    {
      sqrt = (Ne*pow(3E8, 2)*Ee)/(2*M_PI*R*pow(eps, 0.5)*pow(freq, 2.5)*(k*R));
    }
  else
    {
      sqrt = (Ne*pow(3E8, 2)*Ee)/(2*M_PI*R*pow(eps, 0.5)*pow(freq, 2.5)*(1.00-pow(M_E, -k*R)));
    }
  //ans = pow(sqrt, 0.5);
  return sqrt;
}

/*Fn to get max allowed dx from power constraints */
double maxdx_p_new2(double Ne, double Ee, double eps, double freq, double k, double R, double E_low)
{
  double sqrt, ans, param;
  param = (1.000-pow(M_E, -k*R));

  if (param == 0.0)
    {
      sqrt = (Ne*pow(3.0E8, 2.0)*(Ee-E_low))/(2*M_PI*R*pow(eps, 0.5)*pow(freq, 2.5)*(k*R));
    }
  else
    {
      sqrt = (Ne*pow(3.0E8, 3.0)*(Ee-E_low))/(2*M_PI*R*pow(eps, 0.5)*pow(freq, 2.5)*(1.00-pow(M_E, -k*R)));
    }
  ans = pow(sqrt, 0.5);                                                                                              
  return ans;//sqrt;
}



/*Fn to get max allowed dx from power constraints */
double maxdx_pow(double Ne, double Ee, double eps, double freq, double k, double R, double dx)
{
  double sqrt, param;
  param = (1.000-pow(M_E, -k*R));

  if (param == 0.0)
    {
      sqrt = (Ne*pow(3E8, 3)*Ee)/(2*M_PI*R*R*pow(eps, 0.5)*pow(freq, 2.5)*(k*dx));
    }
  else
    {
      sqrt = (Ne*pow(3E8, 3)*Ee)/(2*M_PI*R*R*pow(eps, 0.5)*pow(freq, 2.5)*(1.00-pow(M_E, -k*dx)));
    }
  //ans = pow(sqrt, 0.5);                                                                                              
  return sqrt;
}



/* fn basing the maximum dx on the larmor eqn */
double dx_P_guess(double E_e, double B, double beta, double gamma)
{
  double ans;
  ans = (3*4*M_PI*1E-7*E_e)/(2*6.65E-29*B*B*beta*beta*gamma*gamma);
  return ans;
}


/* fn to compute dx based on both sync & IC */
double dx_synIC(double P_IC, double P_sync, double beta_bulk, double Ebin_min)
// note: P_IC and P_sync need to be per unit length
{
  double sqrt;
  sqrt = pow((beta_bulk*3.0E8*Ebin_min*1.6E-19/(P_IC+P_sync)), 0.50);
  return sqrt;
}


/* fn basing the maximum dx on the larmor eqn */
double dx_P_guess2(double E_e, double B,double beta, double gamma, double E_e2)
{
  double ans;
  ans = (3.0*4.0*M_PI*1E-7*(E_e-E_e2))/(2.0*6.65E-29*B*B*beta*beta*gamma*gamma);
  return ans;
}

/* fn basing the maximum dx on the larmor eqn and IC loss eqn*/
double dx_SIC(double E_e, double Ee2, double B, double beta, double gamma, double U_gam)
{
  double ans;
  ans = (E_e-Ee2)*(3.0/(4.0*6.65E-29*beta*beta*gamma*gamma))*(((2.0*4.0*M_PI*1.0E-7)/(B*B)) + (1.0/U_gam));
  return ans;
}


/* Fn to find a dx using a simple formula to simplify main code */
double dx_P_simple(double Ne, double E_e, double P_tot)
{
  double ans;
  ans = (Ne*E_e*3E8)/P_tot;
  return ans;
}


/* Below are the functions used to compute the IC component of the SED */
double h=6.63E-34, me=9.11E-31, c=3E8, r_c=3.86E-13;


/*fn to model to 4momentum in eqn 2.32c */
double fourmom_photon(double E_ph, double phi2) //input in eV please
{
  double ans;
  ans = 1.0/(1.0+(E_ph/(0.511E6))*(1.0-cos(phi2))); //paper has a plus sign, should it be minus?
  return ans;
}


/*fn to  find the solid angle (2.33) and (2.35) */
double dsolid_angle(double angle, double dangle) //can use for theta or phi2
{
  double ans;
  ans = M_PI*fabs(sin(angle))*dangle;
  return ans;
}


double KN_crosssec(double E_ph_p, double phi2)
{
  double ans;
  ans = 0.5*(1.0/137.0)*(1.0/137.0)*3.86E-13*3.86E-13*pow(fourmom_photon(E_ph_p, phi2), 2.0)*(fourmom_photon(E_ph_p, phi2)+1.0/(fourmom_photon(E_ph_p, phi2)) -1.0 + pow(cos(phi2), 2.0));
  return ans;
}


/* KN cross sec pg237 eqn 9.28 in Longair */
double KN_longair(double E) //E needs to be in eV
{
  double x, re, ans;
  re = 2.81E-15;
  x = E/0.511E6;
  ans = M_PI*re*re*(1/x)*((1.0-2*(x+1.0)/(x*x))*log(2.0*x+1.0) + 0.5 +(4.0/x)-1.0/(2.0*(2.0*x+1.0)*(2.0*x+1.0)));
  return ans;
}


double phden_thin(double L, double freq, double R, double dx)
{
  double ans;
  ans = L/(2.0*M_PI*6.63E-34*freq*R*dx);
  return ans;
} 

/* Fn to get photon density for optically thick kR >1 case */
double phden_thick(double freq, double E_e)
{ 
  double ans, param;
  param = (pow(M_E, (6.63E-34*freq/E_e))-1.0);
  if (param == 0.0)
    {
      ans = (8.0*M_PI*freq*freq*E_e)/(pow(3.0E8, 3.0)*(6.63E-34*freq));
    }
  else
    {
      ans = (8.0*M_PI*freq*freq)/(pow(3.0E8, 3.0)*(pow(M_E, (6.63E-34*freq/E_e))-1.0));
    }
  return ans;
}



/* Fns to explicity do the 4-vector maths for the IC code */

/* Lorentz transform photon energy*/
double pluLT(double E, double gamma, double beta,  double theta)
{
  double ans;
  ans = E*(gamma*(1.0+beta*cos(theta)));
  return ans;
}

/*Inverse Lorentz transform photon energy */
double minLT(double E, double gamma, double beta, double theta)
{
  double ans;
  ans = E*gamma*(1.0-beta*cos(theta));
  return ans;
}


/* transform angles into rest frame */
double angLT(double theta, double gamma, double beta)
{
  double ans;
  ans = asin((sin(theta)/(gamma*(1.0+beta*cos(theta)))));
  return ans;
}


/* inverse transform angles: need two thetas to get back to original frame */
double invangLT(double theta, double gamma, double beta)
{
  double ans;
  ans = asin((sin(theta))/(gamma*(1.0-beta*cos(theta))));
  return ans;
}


/* Fn to compute the photon energy change (loss) in the scattering frame */
double elecscatt(double E_i, double angle)
{
  double ans;
  ans = E_i/(1.0+(E_i/0.511E6)*(1.0-cos(angle)));
  return ans;
}


/* Fn to define the weight function for IC interactions */
double weight_fn(double Ne, double KN, double n_ph, double beta, double theta, double phi2)
{
  double ans;
  ans = (Ne/3.0E8)*KN*(n_ph/(4.0*M_PI))*3.0E8*(1.0+beta*cos(theta))*M_PI*fabs(sin(phi2))*M_PI*fabs(sin(theta));
  return ans;
}



/*Fn to directly compute the upscattering IC fraction assuming some angles */
double IC_onestep(double beta, double gamma, double Eph,  double theta1, double theta2, double alpha)
{
  //theta1 is the angle between elec vel and ph vel initially
  //theta2 between Ee and E_ph'
  //alpha between E_ph and E_ph'
  double ans;
  ans = (1.0-beta*cos(theta1))/((1.0-beta*cos(theta2)+(Eph/(gamma*0.511E6))*(1.0-cos(alpha))));
  return ans;
}



/*Fn to compute the IC losses */
double IC_loss(double beta, double gamma, double Urad)
{
  double ans;
  ans = (4.0/3.0)*6.65E-29*3E8*beta*beta*gamma*gamma*Urad;
  return ans;
}


/*Fn to compute the total set power within one light second */
double P_jet(double B, double gamma, double R, double beta)
{
  double Pow;
  Pow = (B*B/(2.0*4.0*M_PI*1E-7))*3.0E8*M_PI*R*R*gamma*gamma*(1.0+(beta*beta/3.0));
  return Pow;
}


/*Fn to calculate EFFECTIVE B-field due to doppler depolarisation */
double DD_Beffective(double a, double b, double c, double v1, double v2, double v3, double Gamma)
{

  double Blength = sqrt(a*a+b*b+c*c);
  double vlength = sqrt(v1*v1+v2*v2+v3*v3);
  a = a / (Blength); //Normalising B-field vectors
  b = b / (Blength);
  c = c / (Blength);
  v1 = sqrt(1 - 1/(Gamma*Gamma))* v1 / (vlength);
  v2 = sqrt(1 - 1/(Gamma*Gamma))* v2 / (vlength);
  v3 = sqrt(1 - 1/(Gamma*Gamma))* v3 / (vlength);

  double q1 = a + v1*c-v3*a - (Gamma/(1+Gamma))*(v1*a+v2*b+v3*c)*v1;
  double q2 = b + v2*c - b*v3 - (Gamma/(1+Gamma))*(v1*a+v2*b+v3*c)*v2;
  //double q3 = c - (Gamma/(1+Gamma))*(v1*a+v2*b+v3*c)*v3;
  double e1 = -q2 / sqrt(q1*q1+q2*q2);
  double e2 =  q1 / sqrt(q1*q1+q2*q2); //to get vector perp to this (the effective B), swap x and y components and negate new y
  //B_eff = [e2,-e1]
  double Proj_theta_Beff = atan2(-e1,e2);


  return Proj_theta_Beff;

}

/*Fn to calculate EFFECTIVE 3B-field due to doppler depolarisation, now also including z component along our line of sight for IC use */
/* takes B-field and rotates it along plane containing velocity and line of sight by angle defined through Gamma and angle between velocity and line of sight*/
void DD_3Beffective(double a, double b, double c, double v1, double v2, double v3, 
                    double Gamma, double *B_effective)
{

  double Blength = sqrt(a*a+b*b+c*c);
  if (Blength == 0){
    *B_effective = 0;
    *(B_effective+1) = 0;
    *(B_effective+2) = 0;
    return;
  }

  double vlength = sqrt(v1*v1+v2*v2+v3*v3);
  double beta = sqrt(1-1/(Gamma*Gamma));
  double B1, B2, B3, U, V;
  a = a / (Blength); //Normalising B-field vectors
  b = b / (Blength);
  c = c / (Blength);
  v1 = v1 / (vlength);
  v2 = v2 / (vlength);
  v3 = v3 / (vlength);


  double th = acos((v3-beta)/(1-beta*v3)) - acos(v3);
  U = -v2/sqrt(v1*v1+v2*v2);
  V = v1/sqrt(v1*v1+v2*v2);
  B1 = U*U*a + V*V*cos(th)*a + U*V*(1-cos(th))*b + V*sin(th)*c;
  B2 = U*V*(1-cos(th))*a + V*V*b + U*U*cos(th)*b - U*sin(th)*c;
  B3 = -V*sin(th)*a + U*sin(th)*b + (U*U+V*V)*cos(th)*c;

  Blength = sqrt(B1*B1+B2*B2+B3*B3);

  *B_effective = B1/Blength;
  *(B_effective+1) = B2/Blength;
  *(B_effective+2) = B3/Blength;

}

int rand_lim(int limit) {
  /* return a random number between 0 and limit inclusive.
                                                                */

  int divisor = RAND_MAX/(limit+1);
  int retval;

  do {
      retval = rand() / divisor;
  } while (retval > limit);

  return retval;
  }

double theta_X(double theta_tot, double theta_circ) {
  /*    */
  double theta_x = asin(sin(theta_tot)/sqrt(1+1/(tan(theta_circ)*tan(theta_circ))));
  return theta_x;
  }

double theta_Y(double theta_tot, double theta_circ) {
  /*    */
  double theta_y = asin(sin(theta_tot)*cos(theta_circ));
  return theta_y;
  }

double Z_perpperp(double cosk, double phik, double Btheta) {
  /* collision factor Z for Compton emission perpendicular to B-field with target photons polarised perp to B _field  */
  return pow((sin(Btheta)*cosk-cos(Btheta)*sqrt(1-cosk*cosk)*cos(phik)-(1-cosk*cosk)*sin(phik)*sin(phik)*sin(Btheta)/(1-cosk))/sqrt(pow(sqrt(1-cosk*cosk)*sin(phik),2) + pow(cos(phik)*cos(Btheta)*sqrt(1-cosk*cosk) - cosk*sin(Btheta),2)),2);
  }

double Z_perppara(double cosk, double phik, double Btheta) { //still to do
  /* collision factor Z for Compton emission perpendicular to B-field with target photons polarised perp to B _field  */
  double e_mag = sqrt(pow((1-cosk*cosk)*pow(cos(phik),2)*sin(Btheta) + sqrt(1-cosk*cosk)*cos(phik)*cos(Btheta)*cosk-sin(Btheta),2) +
                pow((1-cosk*cosk)*cos(phik)*sin(phik)*sin(Btheta) + sqrt(1-cosk*cosk)*sin(phik)*cos(Btheta)*cosk,2) +
                pow(cosk*sqrt(1-cosk*cosk)*cos(phik)*sin(Btheta) + (cosk*cosk - 1)*cos(Btheta),2) );
  return pow(((1-cosk*cosk)*cos(phik)*sin(phik)*sin(Btheta) + sqrt(1-cosk*cosk)*sin(phik)*cos(Btheta)*cosk + sqrt(1-cosk*cosk)*sin(phik)*(cosk*sqrt(1-cosk*cosk)*cos(phik)*sin(Btheta) + (cosk*cosk - 1)*cos(Btheta))/(1-cosk))/e_mag,2);
  }

double Z_paraperp(double cosk, double phik, double Btheta) {
  /* collision factor Z for Compton emission perpendicular to B-field with target photons polarised perp to B _field  */
  return pow((sqrt(1-cosk*cosk)*sin(phik)*cos(Btheta)-(1-cosk*cosk)*sin(phik)*cos(phik)*sin(Btheta)/(1-cosk))/sqrt(pow(sqrt(1-cosk*cosk)*sin(phik),2) + pow(cos(phik)*cos(Btheta)*sqrt(1-cosk*cosk) - cosk*sin(Btheta),2)),2);
  }

double Z_parapara(double cosk, double phik, double Btheta) {//still to do
  /* collision factor Z for Compton emission perpendicular to B-field with target photons polarised perp to B _field  */
  double e_mag = sqrt(pow((1-cosk*cosk)*pow(cos(phik),2)*sin(Btheta) + sqrt(1-cosk*cosk)*cos(phik)*cos(Btheta)*cosk-sin(Btheta),2) +
                pow((1-cosk*cosk)*cos(phik)*sin(phik)*sin(Btheta) + sqrt(1-cosk*cosk)*sin(phik)*cos(Btheta)*cosk,2) +
                pow(cosk*sqrt(1-cosk*cosk)*cos(phik)*sin(Btheta) + (cosk*cosk - 1)*cos(Btheta),2) );
  return pow(((1-cosk*cosk)*pow(cos(phik),2)*sin(Btheta) + sqrt(1-cosk*cosk)*cos(phik)*cos(Btheta)*cosk-sin(Btheta) + sqrt(1-cosk*cosk)*cos(phik)*(cosk*sqrt(1-cosk*cosk)*cos(phik)*sin(Btheta) + (cosk*cosk - 1)*cos(Btheta))/(1-cosk))/e_mag,2);
}

//just have one function for Z which takes inputs to choose polarisation vectors of incoming and outgoing:
double Z_e(double *_e, double *e, double *k, int perp) {

  double e_length = sqrt(e[0]*e[0]+e[1]*e[1]+ e[2]*e[2]);
  e[0]=e[0]/e_length, e[1]=e[1]/e_length, e[2]=e[2]/e_length;
  double k_length = sqrt(k[0]*k[0]+k[1]*k[1]+ k[2]*k[2]);
  k[0]=k[0]/k_length, k[1]=k[1]/k_length, k[2]=k[2]/k_length;

  if (perp==1) {
    _e[0] = e[0] + e[2]*k[0]/(1-k[2]), _e[1] = e[1] + e[2]*k[1]/(1-k[2]);
    //_e[0] = 0, _e[1] = 1; //original fixed setting
  } else {
    _e[0] = -e[1] - e[2]*k[1]/(1-k[2]), _e[1] = e[0] + e[2]*k[0]/(1-k[2]);
    //_e[0] = -1, _e[1] = 0; //original fixed setting
  }
  double _e_length = sqrt(_e[0]*_e[0]+_e[1]*_e[1]);
  _e[0]=_e[0]/_e_length, _e[1]=_e[1]/_e_length;

  return pow(e[0]*_e[0] + e[1]*_e[1] + (k[0]*_e[0] + k[1]*_e[1])*(e[2])/(1-k[2]),2);
}

double Sigma_1(double E_elecs,double F_min, double dN_dE, double dEe ) {
  /* Sigma1 energy integral */
  return (dN_dE*dEe*1.6E-19*F_min*(pow(F_min/(E_elecs*1.6E-19),2) - pow(E_elecs*1.6E-19/F_min,2) + 2)/pow(E_elecs*1.6E-19,4));
  }

double Sigma_2(double E_elecs,double F_min, double dN_dE, double dEe ) {
  /* Sigma1 energy integral */
  return (dN_dE*dEe*1.6E-19*pow(pow(E_elecs*1.6E-19,2) - F_min*F_min,2)/(pow(E_elecs*1.6E-19,6)*F_min));
  }


/* Prepends t into s. Assumes s has enough space allocated
 * ** for the combined string.
 * */
void prepend(char* s, const char* t)
{
    size_t len = strlen(t);
    size_t i;

    memmove(s + len, s, strlen(s) + 1);

    for (i = 0; i < len; ++i)
    {
        s[i] = t[i];
    }
}
