#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "jet_fns.h"

/* Program to store some useful Lorentz transform and jets fns */

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


/* Fn to return a logarithmically spaced array */
void log10space(float array[], float start, float end)
{
  float endstart;
  int no;
  endstart = end/start;
  //  printf("logendstart: %.2f \n", log10(endstart));
  no = log10(endstart);
  //printf("The no of pts is: %d \n", no);

  float element;//array[no],  element;
  int j;
  element = start;
  //printf("no %.2d", no);
  for (j=0; j<=no; j+=1)//element<=end; j++)
    {
      array[j]=element;
      element *= 10;
      //  printf("\n loop no is: %d\t%.f", j,element);
 
   }
  //  return 0;//(array);
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
      if (array[i]<min && array[i]!=0.0) //needs to be non zero or code won't work!
	{
	  //printf("%.2e\t%s\t%.2e\t%d \n", array[i], "<", min, element);
	  min = array[i];
	  element = i;
	}
    }
  return element;
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