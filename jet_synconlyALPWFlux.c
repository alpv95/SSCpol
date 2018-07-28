#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "jet_fns.h"
#include <time.h>

/* Program model synchrotron emission with user friendly interface

last edit: 19/10/16
created: 19/10/16

outputs: correcty obtains synchrotron bump

IC status:
-energy upscatter appears to work as no upscatter > 4*gamma^2 possible

*/

//define some physical constants
double me_eV=0.511E6; //electron rest mass



//define some looping parameters
int array_size=50; // sets the number of synchrotron & IC pts
double array_size_d=50.0; // use to set log ratios, should be the same as above line but with .0
int i, l, m, n, o, p; //some looping parameters
double w, x, y, z, a, b, c, d, q;
int dx_set; // use to define smallest non zero population
int counter;

int main(int argc,char* argv[]) //argc is integer number of arguments passed, argv[0] is program name, argv[1..n] are arguments passed in string format
{
    //printf("%s",argv[1]); //wanna figure out how to get inputs to come in not as strings
    FILE *OW_jfile, *test, *FGfile;
    OW_jfile = fopen("OW_jfile.txt","r");
    test = fopen("test.txt","r");
    FGfile = fopen("FG.txt","r");

    double FG[75][3];
    double W_jarray[100];
    /*double L_jarray[30];
    double Barray[30];
    double E_maxarray[30];
    double theta_obsarray[30];
    double alphaarray[30];
    double gamma_bulkarrayMax[30];
    double gamma_bulkarrayMin[30];*/
    
    fscanf(test, "%d", &p);
    counter = p;
    
    for (m=0; m<100; m++)
    {
        fscanf(OW_jfile, "%lf\n", &w);
        W_jarray[m] = w;
        if (m<75)
        {
           fscanf(FGfile, "%lf\t%lf\t%lf\n", &q, &c, &d);
           FG[m][0] = q;
           FG[m][1] = c;
           FG[m][2] = d;
        }
        
    }
    //enter the parameters from the jets paper
    double W_j=2E37;  //3e37//W_jarray[counter]; //W jet power in lab frame 7.63E36. Should be OBSERVED POWER
     double L_jet = 5E20;// 6E11;//1E19;//6E12;// 1E19;// 6E12;//6E12; //1.0E19; //length in m in the fluid frame
     double B=3E-3, B0=3E-3;  //8e-5//0.000268; //B-field at jet base
     double R0 = 0.0, R=0.0; //4.532E13;//7.32E13; // Radius of the jet at the base 3.32 works fairly well 
    double R_prev=0.0, B_prev=0.0; //changing parameters of the jet-initialise. R prev corrects for increasing jet volume 
    double E_min = 5.11E6; // Minimum electron energy 
    double E_max=32E9;//50e9//8.1E9//5.0E9;//5.60E9; // Energy of the ECO in eV 
    double alpha=1.95;//1.9//1.95//2.000001; // PL index of electrons
     double theta_open_p = 50.0;//50//*(M_PI/180.0); // opening angle of the jet in the fluid frame 
    double theta_obs; //4//3.0;//*(M_PI/180.0); // observers angle to jet axis in rad 
    sscanf(argv[4], "%lf", &theta_obs);
    double gamma_bulk=17.5;//pow(10,(log10(W_j)*0.246-8.18765 + 0.09)); //final additive constant to make sure highest is 40 and lowest is 5//12.0; // bulk Lorentz factor of jet material
    printf("%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n", W_j, gamma_bulk, theta_obs, B, B0, E_max, alpha);



  //define radius at jet base
  R0 = R_0(W_j*(3.0/4.0), 1.0, gamma_bulk, B0);//7.32E13;
  R = R_0(W_j*(3.0/4.0), 1.0, gamma_bulk, B0);//7.32E13;
  printf("Radius at jet base: %.5e \n", R);

    //define some files to store output data
  FILE *Nfile, *freqrange, *freqfile, *opfile, *powfile, *basicdata, *elecpop, *jfile, *keyparams, *ICfile, *ICfreqfile, *Pperpfile, *Pparafile, *Proj_Bfile, *block_thetafile;
    Nfile = fopen("Ndata.txt", "w"); //store electron populations
    freqrange = fopen("freqrange.txt", "w");//store frequency bin boundaries
    freqfile = fopen("critfreqs.txt", "w");//store critical frequencies
    opfile = fopen("kvalues.txt", "w");//store opacity values for each jet section
    powfile = fopen("secpow.txt", "w");//store emitted synchrotron power (freq fn) for each jet section
    basicdata = fopen("basicdata.txt", "w"); //store B, x, R etc
    elecpop = fopen("elecpop.txt", "w"); //store electron populations at various points down the jet
    jfile = fopen("jfile_re.txt", "w");
    keyparams = fopen("keyparams.txt", "w"); //store Lj, gamma_bulk and theta_obs
    ICfile = fopen("ICfile.txt", "w"); //ALP
    ICfreqfile = fopen("ICfreqfile.txt", "w"); //ALP
    Pperpfile = fopen("Pperpfile.txt", "w"); //polarisations
    Pparafile = fopen("Pparafile.txt", "w");
    Proj_Bfile = fopen("Proj_Bfile.txt", "w"); //projected B field onto plane of the sky for each section
    block_thetafile = fopen("block_thetafile.txt","w"); //saves the angle to the line of sight of each block
    fprintf(keyparams, "\t%.5e\t%.5e\t%.5e \n", L_jet, gamma_bulk, theta_obs);

    //define some useful parameters to gauge progress
    int nSteps=0; //how many sections the jet has broken down into
    double x=0; //progress along the jet
    double dx; //incremenet of step length
    double eps; //stores the eqn 2.18b
    double beta_bulk = lortovel(gamma_bulk);
    double doppler_factor = doppler(beta_bulk, deg2rad(theta_obs));
    int intermed_param=0; //allows intermediate population to be determined to compare to paper

    //define some parameters for determining the dx of each jet section
    double dx_R; //dx where R_new = 1.05*R_old
    double dx_P; //dx to ensure smallest bin population not depleted

    //some parameters used to define the population of electrons
    double A_elecs[array_size]; //stores PL coefficient
    //double E_elecs[array_size]; // stores the energy in eV of each electron
    double Ne_e[array_size]; // no of electrons in each bin
    double Ne_orig[array_size]; // at jet base: allows comparison later
    double Ne_intermed[array_size]; // intermediate population
    double dN_dE[array_size]; // no of electrons per J in each bin
    double dN_dE_orig[array_size];
    double dEe[array_size]; //size of bin widths in eV
    double Ne_losses[array_size]; //losses per section
    double Ne_gains[array_size]; //gains per section
    double A_orig; // PL prefactor initially the same for all electrons
    double Ue_jetbase = 0.0; // energy in particles at conical region base
    double UB_jetbase = 0.0; //energy in particles at the jet base

    // parameters used to initialise electron population
    //double Ee_max=4.0E12; //sets the maximum electron energy in eV
    //double E0_mid; //defines the mid bin value of the first bin-need lower error bar equal to m_e in eV
    //double ratio=1.32; //logarithmic ratio of electron bins-done by hand
    //E0_mid = 0.511E6*(pow(ratio, 0.5)); //sets lower error as electron rest mass
    //elecErange(E_elecs, E0_mid, Ee_max, array_size_d); // defines bin midpts from E0mid to Ee_max in array E_elecs
    //ratio = pow((Ee_max/E0_mid), (1/array_size_d)); //redefines the exact ratio
    //double E_elec_min[array_size], E_elec_max[array_size]; //these map onto critical frequencies

    // parameters used to initialise electron population
    double Ee_max=4.0E12; //sets the maximum electron energy in eV
    double ratio= pow((Ee_max/E_min), (1/array_size_d)); //logarithmic ratio of electron bins: USED FOR BOUNDS
    double E_bounds[array_size+1], indices[array_size+1];
    arange(indices, array_size+1); //as np.arange
    
    
    for (i=0; i<array_size+1; i++)
    {
    E_bounds[i] = E_min*(pow(ratio, indices[i]));
    //printf("indices, bounds: %.5e\t%.5e \n", indices[i], E_bounds[i]); 
    }
    
    //now sort the bounds into max, min and E_elecs
    double E_elecs[array_size]; // stores the energy in eV of each electron
    double E_elec_min[array_size], E_elec_max[array_size]; //these map onto critical frequencies, should be in eV
    
    for (i=0; i<array_size; i++)
    {
    E_elec_min[i] = E_bounds[i];
    E_elec_max[i] = E_bounds[i+1];
    E_elecs[i] = pow(E_elec_min[i]*E_elec_max[i], 0.5); 
    }

    //parameters based on individual electron bins
    double gamma_e[array_size]; // store Lorentz factors
    double beta_e[array_size]; // v/c

    //parameters used in determining synchrotron emission
    double f_c[array_size]; //store the critical frequencies
    double f_c_IC[array_size]; //ALP
    double dfreqs[array_size];// store the bin widths in critical frequency
    double f_em_min[array_size], f_em_max[array_size], f_em_min_IC[array_size], f_em_max_IC[array_size]; // store bin boundaries for critical frequencies
    double j[array_size]; // emissivity per Hz per unit volume
    double k[array_size]; // opacity, units of m^-1
    double Ps_per_m[array_size]; // Synchrotron power per unit m
    double Ps_per_m_elec[array_size];
    double Ps_per_m_test[array_size];
    double Sync_losses[array_size];// store total electron energy losses
    double P_single[array_size];// power emitted by a single electron of energy E

    //Fill electron arrays with initial parameters
    for (i=0; i<array_size; i++)
    {
        E_elec_min[i] = E_elecs[i]*(1.0/pow(ratio, 0.5)); //lower electron eV bounds
        E_elec_max[i] = E_elecs[i]*pow(ratio, 0.5); //upper electron eV bounds
        dEe[i] = E_elec_max[i]-E_elec_min[i]; // sets bin widths in eV
        A_elecs[i] = A_PL(alpha, W_j, gamma_bulk, E_min*1.6E-19, E_max*1.6E-19, 1.0); //1 assumes equipartition. DIVIDED W_j by gamma^2 for jet frame
        dN_dE[i] = electron_PL(A_elecs[i], alpha, E_elecs[i]*1.6E-19, E_max*1.6E-19);
        dN_dE_orig[i] = electron_PL(A_elecs[i], alpha, E_elecs[i]*1.6E-19, E_max*1.6E-19);
        Ne_e[i] = electron_PL(A_elecs[i], alpha, E_elecs[i]*1.6E-19, E_max*1.6E-19)*dEe[i]*1.6E-19;
        Ne_orig[i] = electron_PL(A_elecs[i], alpha, E_elecs[i]*1.6E-19, E_max*1.6E-19)*dEe[i]*1.6E-19; //allows comparison of final to initial electron population

        //now calculate the initial critical frequencies
        gamma_e[i] = E_elec_min[i]/me_eV; //works as all in eV
        f_em_min[i] = f_crit(gamma_e[i], B); //lower bin values for f_crits
        f_em_min_IC[i] = gamma_e[i] * gamma_e[i] * f_em_min[i]; //* (6.63E-34 /(1.6E-19))

        gamma_e[i] = E_elec_max[i]/me_eV; //redefine for upper boundary
        f_em_max[i] = f_crit(gamma_e[i], B);
        f_em_max_IC[i] = gamma_e[i] * gamma_e[i] * f_em_max[i]; //* (6.63E-34 /(1.6E-19))
        
        gamma_e[i] = E_elecs[i]/me_eV; //mid points
        beta_e[i] = lortovel(gamma_e[i]); // beta taken at midpts

        f_c[i] = f_crit(gamma_e[i], B); //assume all radiation occurs at the critical frequency
        f_c_IC[i] =gamma_e[i] * gamma_e[i] * f_c[i]; //ALP
        dfreqs[i] = f_em_max[i]-f_em_min[i]; //obtain frequency bin widths
        
        Ue_jetbase += Ne_e[i]*E_elecs[i]*1.6E-19; 

        fprintf(freqrange,  "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e \n", f_em_min[i], f_em_max[i], f_c[i], dfreqs[i], f_em_min_IC[i], f_em_max_IC[i], f_c_IC[i]);//saves the frequency bins!

    }//ends electron population definition.
    A_orig = A_elecs[0]; //all elements the same to begin with


    //print the energies to test equipartition
    UB_jetbase = (4.0/3.0)*M_PI*R0*R0*3.0E8*B0*B0/(2*4*M_PI*1E-7);
    //printf("Obs power, particles, B-fields: %.5e\t%.5e\t%.5e\t%.5e \n", W_j/(gamma_bulk*gamma_bulk), Ue_jetbase, UB_jetbase, Ue_jetbase+UB_jetbase);
    
    //redefine electron population to ensure equipartition
    for (i=0; i<array_size; i++)
    {

        A_elecs[i] = A_PL(alpha, W_j, gamma_bulk, E_min*1.6E-19, E_max*1.6E-19, 1.0)*UB_jetbase/Ue_jetbase; //1 assumes equipartition. DIVIDED W_j by gamma^2 for jet frame
        dN_dE[i] = electron_PL(A_elecs[i], alpha, E_elecs[i]*1.6E-19, E_max*1.6E-19);
        Ne_e[i] = electron_PL(A_elecs[i], alpha, E_elecs[i]*1.6E-19, E_max*1.6E-19)*dEe[i]*1.6E-19;
        Ne_orig[i] = electron_PL(A_elecs[i], alpha, E_elecs[i]*1.6E-19, E_max*1.6E-19)*dEe[i]*1.6E-19; //allows comparison of final to initial electron population

    }//ends electron population re-definition.
    
    //recalaculate equipartition fraction
    Ue_jetbase = 0.0; 
    for (i=0; i<array_size; i++)
    {
        Ue_jetbase += Ne_e[i]*E_elecs[i]*1.6E-19; 
    }

    //fudge factors: set bins with fewer than ten electrons to zero: prevents resolution from being limited
    cleanpop(Ne_e, array_size, 10.0);
    cleanpop(dN_dE, array_size, 10.0);
    

    for (i=0; i<array_size; i++)
    {
        A_elecs[i]*=(Ne_e[i]/Ne_orig[i]); //sets A to zero for tiny populations
    }



    //define photon energy density
    double u_rad=0.0;
    double urad_array[array_size]; //convert synchrotron power into radiation field
    
    
    //11/10/16 define some parameters to determine thick/thin jet sections
    double R_eff[array_size]; //effective radius down to which can be seen -> smaller than jet radius if optically thick

    //ALP initialisation
    double Ps_per_m_IC[array_size];
    double Ne_etot;
    double n_rad;
    double IC_losses[array_size];

    //****************Initialise Polarisation variables *****************************************************************************//
    //first set frequency bins for polarised powers to drop into:
    double f_pol[array_size];
    double P_perp[array_size];
    double P_para[array_size];
    double dfreqs_pol[array_size];
    //initialise the variables with appropriate values:
    for (i=0; i<array_size; i++)
    {
        f_pol[i] = f_c[i];
        dfreqs_pol[i] = dfreqs[i];
    }
    //HELICAL B-field: This variable will save the projection of the B-field on the sky for each section of the jet
    //now it is an array to include all the random blocks of B-field
    int n_blocks = 127; //for the marscher model, can have 1,7,19,37,61,91,127 blocks
    double Proj_theta_B[n_blocks];
    //defining how far along the helix we are
    int helix_counter;
    sscanf(argv[2], "%d", &helix_counter); //taking input from command line how many 'days' we have been observing
    //parameters for helix
    double t_helix = helix_counter*M_PI/16; //helix parameter x=rcos(t+phi) y=rsin(t+phi) z =ct (r will be radius of jet)
    //Pi constant decided how quickly we sample along the helix -> timescale
    double phi_helix = M_PI/3; //phase shift, just an initial condition
    double c_helix = R0; //speed (ie number of coils per unit distance) -- high c => spaced out coils, low c => densely packed coils
    //choosing random vectors (a_0,b_0,c_0) in unit sphere for random blocks initial B directions:
    srand((unsigned)time(NULL)); //seeds random number generator from computer clock
    //srand(helix_counter); //for keeping seeds the same
    double a_0[n_blocks-1];
    double b_0[n_blocks-1];
    double c_0[n_blocks-1];
    double B_length;
    for (i=0; i<(n_blocks); i++){ //these are the random B-field vectors in each block as if we are looking straight down jet (have to rotate by theta_obs)
      a_0[i] = rand_lim(200);
      b_0[i] = rand_lim(200);
      c_0[i] = rand_lim(200);
      a_0[i]=(a_0[i]-100)/100;
      b_0[i]=(b_0[i]-100)/100;
      c_0[i]=(c_0[i]-100)/100;
      B_length = sqrt(a_0[i]*a_0[i]+b_0[i]*b_0[i]+c_0[i]*c_0[i]);
      a_0[i] = a_0[i] / (B_length); //Normalising B-field vectors
      b_0[i] = b_0[i] / (B_length);
      c_0[i] = c_0[i] / (B_length);
    }
    int EVPA_rotation;
    sscanf(argv[1], "%d", &EVPA_rotation); //taking input from command line if true begin EVPA rotation, if false do not
    int DopDep;
    sscanf(argv[3], "%d", &DopDep); //taking input from command line if true activate DopDep, if false do not

    //now the Bfield vectors are rotated about both the y and x axis each block with a different theta_obs
    //these give the angles each block is at with respect to the direction of the jet:
    //theta_x => rotation angle about x axis, theta_y => rotation angle about y axis
    double th = 2*atan(tan(deg2rad(theta_open_p))/gamma_bulk); //full jet diameter opening angle in lab frame
    /* double theta_x[49] = {0,0,0,0,0,                            // 5 5 5 configuration for central helical field
                     th/9,th/9,th/9,th/9,th/9,
                     -th/9,-th/9,-th/9,-th/9,-th/9,-th/9,-th/9,0,0,th/9,th/9,
                     2*th/9,2*th/9,2*th/9,2*th/9,2*th/9,2*th/9,
                     3*th/9,3*th/9,3*th/9,3*th/9,3*th/9,
                     4*th/9,4*th/9,4*th/9,
                     -2*th/9,-2*th/9,-2*th/9,-2*th/9,-2*th/9,-2*th/9,
                     -3*th/9,-3*th/9,-3*th/9,-3*th/9,-3*th/9,
                     -4*th/9,-4*th/9,-4*th/9};
    //7+7+7+6+6+5+5+3+3 configuration
    double theta_y[49] = {0,th/7,2*th/7,-th/7,-2*th/7,
                      0,th/7,2*th/7,-th/7,-2*th/7,
                      0,th/7,2*th/7,-th/7,-2*th/7,-3*th/7,3*th/7,3*th/7,-3*th/7,3*th/7,-3*th/7,
                      sqrt(th*th-2*(2*th/9)*(2*th/9))/6,2*sqrt(th*th-2*(2*th/9)*(2*th/9))/6,3*sqrt(th*th-2*(2*th/9)*(2*th/9))/6,-sqrt(th*th-2*(2*th/9)*(2*th/9))/6,-2*sqrt(th*th-2*(2*th/9)*(2*th/9))/6,-3*sqrt(th*th-2*(2*th/9)*(2*th/9))/6,
                      0,sqrt(th*th-2*(3*th/9)*(3*th/9))/5,2*sqrt(th*th-2*(3*th/9)*(3*th/9))/5,-sqrt(th*th-2*(3*th/9)*(3*th/9))/5,-2*sqrt(th*th-2*(3*th/9)*(3*th/9))/5,
                      0,sqrt(th*th-2*(4*th/9)*(4*th/9))/3,-sqrt(th*th-2*(4*th/9)*(4*th/9))/3,
                      sqrt(th*th-2*(2*th/9)*(2*th/9))/6,2*sqrt(th*th-2*(2*th/9)*(2*th/9))/6,3*sqrt(th*th-2*(2*th/9)*(2*th/9))/6,-sqrt(th*th-2*(2*th/9)*(2*th/9))/6,-2*sqrt(th*th-2*(2*th/9)*(2*th/9))/6,-3*sqrt(th*th-2*(2*th/9)*(2*th/9))/6,
                      0,sqrt(th*th-2*(3*th/9)*(3*th/9))/5,2*sqrt(th*th-2*(3*th/9)*(3*th/9))/5,-sqrt(th*th-2*(3*th/9)*(3*th/9))/5,-2*sqrt(th*th-2*(3*th/9)*(3*th/9))/5,
                      0,sqrt(th*th-2*(4*th/9)*(4*th/9))/3,-sqrt(th*th-2*(4*th/9)*(4*th/9))/3}; */

    //number of rings = 4 (not including middle), number of cells = 61
    int n_rings = 6; //up to 6 rings possible atm
    double theta_r[n_blocks];                // Marscher configuration, rings of turbulent cells, middle 3 rings turbulent
    double theta_phi[n_blocks];
    //Marscher config, calculates theta_r and theta_phi for each block in the concentric circles
    for (i=0; i<n_blocks; i++){
        if (i<1){
            theta_r[i] = 0;
            theta_phi[i] = 0;
        }
        else if (i>0 && i<7){
            for (l=0; l<6; l++){
                theta_r[i] = th/(2*n_rings);
                theta_phi[i] = -M_PI+l*2*M_PI/6;
            }
        }
        else if (i>6 && i<19){
            for (l=0; l<12; l++){
                theta_r[i] = 2*th/(2*n_rings);
                theta_phi[i] = -M_PI+l*2*M_PI/12;
            }
        }
        else if (i>18 && i<37){
            for (l=0; l<18; l++){
                theta_r[i] = 3*th/(2*n_rings);
                theta_phi[i] = -M_PI+l*2*M_PI/18;
            }
        }
        else if (i>36 && i<61){
            for (l=0; l<24; l++){
                theta_r[i] = 4*th/(2*n_rings);
                theta_phi[i] = -M_PI+l*2*M_PI/24;
            }
        }
        else if (i>60 && i<91){
            for (l=0; l<30; l++){
                theta_r[i] = 5*th/(2*n_rings);
                theta_phi[i] = -M_PI+l*2*M_PI/30;
            }
        }
        else if (i>90 && i<127){
            for (l=0; l<36; l++){
                theta_r[i] = 6*th/(2*n_rings);
                theta_phi[i] = -M_PI+l*2*M_PI/36;
            }
        }
    }

    double theta_tot[n_blocks]; //total angle to line of sight of each block (for doppler boosting) in radians
    for (i=0; i<n_blocks; i++){
        theta_tot[i] = acos(cos(theta_r[i])*cos(deg2rad(theta_obs)) + sin(theta_r[i])*sin(deg2rad(theta_obs))*cos(M_PI - fabs(theta_phi[i])));
    }

    for (i=0; i<n_blocks; i++){
        fprintf(block_thetafile, "\t%.5e", theta_tot[i]); //saving the angle to the direction of view of each block to be used in python file for doppler boost
    }
    /*
    double doppler_factor[n_blocks];
    for (i=0; i<n_blocks; i++){
        doppler_factor[i] = doppler(beta_bulk, theta_tot[i]);
    }
    */
    double a_00[n_blocks];
    double b_00[n_blocks];
    double c_00[n_blocks];
    double a_helical[n_blocks]; //simplifies angular shifts for the helical components
    double b_helical[n_blocks];
    double c_helical[n_blocks];


    printf("Beginning jet analysis... \n");
    while (x<L_jet)//for (x=0; x<L_jet; x+=dx)
    {
        //update jet parameters
        eps = epsilon(B);//same for all populations
	//printf("x, R, B; %.5e\t%.5e\t%.5e\n", x, R, B); //matches

        //begin with Synchrotron emission
        for (i=0; i<array_size; i++)
        {
            P_single[i] = larmor(B, beta_e[i], gamma_e[i]); //gives the radiated power by one electron in each bin
            j[i] = j_per_hz(P_single[i], R, dN_dE[i], eps, f_c[i]);//*(R_next/R)*(R_next/R);//gives emissivity per unit frequency

            k[i] = k_new(j[i], eps, f_c[i]); //gives the corresponding opacity
            //Ps_per_m[i] = power_emitted(j[i], k[i], R)*M_PI*R*dfreqs[i]; //power per unit length assuming R roughly constant
            //Ps_per_m[i] = P_single[i]*(2*M_PI*pow(9.11E-31,2)*pow(3E8,2))*dN_dE[i]*dfreqs[i]/(3E8*3*gamma_e[i]*1.6E-19*B);
            //Ps_per_m[i] = P_single[i]*dN_dE[i]*E_elecs[i]*1.6E-19*dfreqs[i] / (2*f_c[i]*3E8);
            //printf("Ps_per_m P_single %.5e\t%.5e\n", Ps_per_m[i],P_single[i]);
            //Sync_losses[i] = Ps_per_m[i];

            //determine whether jet is optically thin or thick -> this part is incomplete
            R_eff[i] = power_emitted(j[i], k[i], R)/j[i]; //depth down to which can be seen in one second


            //Polarisation -> Power emitted parallel and perpendicular to projected B field at each frequency interval

            for (l=0; l<array_size; l++)
            {
                for (m=0; m<75; m++)
                {
                    if ((f_pol[l]/f_c[i]) >= FG[m][0] && (f_pol[l]/f_c[i]) <= FG[m+1][0])
                    {
                        P_perp[l] += (sqrt(3)*M_PI/(16*9.11E-31*pow(3E8,3)))*pow(1.6E-19,3) * B * FG[m+1][1] * dN_dE[i]*E_elecs[i]*1.6E-19 * dfreqs_pol[l] / (1E-7);     //power per unit frequency, have to *dfreq[j]
                        P_para[l] += (sqrt(3)*M_PI/(16*9.11E-31*pow(3E8,3)))*pow(1.6E-19,3) * B * FG[m+1][2] * dN_dE[i]*E_elecs[i]*1.6E-19 * dfreqs_pol[l] / (1E-7) ; //*dx/l_c !!! (this happens further down) and * other constants
                        Ps_per_m_elec[i] += (sqrt(3)*M_PI/(16*9.11E-31*pow(3E8,3)))*pow(1.6E-19,3) * B * (FG[m+1][1]+FG[m+1][2]) * dN_dE[i]*E_elecs[i]*1.6E-19 * dfreqs_pol[l] / (1E-7) ;
                    }
                    else if ((f_pol[l]/f_c[i]) <= FG[0][0] || (f_pol[l]/f_c[i]) >= FG[73][0])
                    {
                        P_perp[l] += 0.0;
                        P_para[l] += 0.0;
                        Ps_per_m_elec[i] += 0.0;
                    }
                }
            }
            Sync_losses[i] = Ps_per_m_elec[i];
            //printf("Ps_per_m Ps_per_me_elec %.5e\t%.5e\n", Ps_per_m[i], Ps_per_m_elec[i]);

        }

        for (l=0; l<array_size; l++){
            Ps_per_m[l] = P_perp[l] + P_para[l];
            //printf("Ps_per_m_test %.5e\n",Ps_per_m_test[l]);
            //photon energy density
            urad_array[l] = 0.0; //initialise to avoid weird C bugs
            urad_array[l] = Ps_per_m[l]/(2.0*M_PI*3.0E8*R);
            u_rad += urad_array[l];
        }



        //----------------- B-field projection onto sky for this section --------------//
        // Now have introduced R dependence on the random blocks as well, works just like for the helical case, just * R_old/R_new on c_0 term
        //now can also include doppler depolarisation, changing the B field to EFFECTIVE B fields which give correct perp E component as if it had
        //been doppler depolarisedx
        //Now each block has its own theta_obs as we are looking inside the conical jet
        for (i=0; i<(n_blocks); i++){
        a_00[i] = a_0[i]*(cos(theta_r[i])+cos(M_PI/2+theta_phi[i])*cos(M_PI/2+theta_phi[i])*(1-cos(theta_r[i]))) + b_0[i]*cos(M_PI/2+theta_phi[i])*sin(M_PI/2+theta_phi[i])*(1-cos(theta_r[i])) + c_0[i]*R0/R*sin(M_PI/2+theta_phi[i])*sin(theta_r[i]);
        b_00[i] = a_0[i]*sin(M_PI/2+theta_phi[i])*cos(M_PI/2+theta_phi[i])*(1-cos(theta_r[i])) + b_0[i]*(cos(theta_r[i])+sin(M_PI/2+theta_phi[i])*sin(M_PI/2+theta_phi[i])*(1-cos(theta_r[i]))) - c_0[i]*R0/R*cos(M_PI/2+theta_phi[i])*sin(theta_r[i]);
        c_00[i] = -a_0[i]*sin(M_PI/2+theta_phi[i])*sin(theta_r[i]) + b_0[i]*cos(M_PI/2+theta_phi[i])*sin(theta_r[i]) + c_0[i]*R0/R*cos(theta_r[i]);
        //helical B-field equivalent has to be inside loop below as it depends on R
        }

        for (i=0; i<(n_blocks); i++){
            //a_helical[i] = -R*sin(t_helix+phi_helix)*cos(theta_y[i]) + sin(theta_y[i])*(R*cos(t_helix+phi_helix)*sin(theta_x[i])+c_helix*sin(theta_x[i]));
            //b_helical[i] = R*cos(t_helix+phi_helix)*cos(theta_x[i]) - c_helix*sin(theta_x[i]);
            //c_helical[i] = R*sin(t_helix+phi_helix)*sin(theta_y[i]) + cos(theta_y[i])*(R*cos(t_helix+phi_helix)*sin(theta_x[i])+c_helix*cos(theta_x[i]));
            a_helical[i] = -R*sin(t_helix+phi_helix)*(cos(theta_r[i])+cos(M_PI/2+theta_phi[i])*cos(M_PI/2+theta_phi[i])*(1-cos(theta_r[i]))) + R*cos(t_helix+phi_helix)*cos(M_PI/2+theta_phi[i])*sin(M_PI/2+theta_phi[i])*(1-cos(theta_r[i])) + c_helix*sin(M_PI/2+theta_phi[i])*sin(theta_r[i]);
            b_helical[i] = -R*sin(t_helix+phi_helix)*sin(M_PI/2+theta_phi[i])*cos(M_PI/2+theta_phi[i])*(1-cos(theta_r[i])) + R*cos(t_helix+phi_helix)*(cos(theta_r[i])+sin(M_PI/2+theta_phi[i])*sin(M_PI/2+theta_phi[i])*(1-cos(theta_r[i]))) - c_helix*cos(M_PI/2+theta_phi[i])*sin(theta_r[i]);
            c_helical[i] = R*sin(t_helix+phi_helix)*sin(M_PI/2+theta_phi[i])*sin(theta_r[i]) + R*cos(t_helix+phi_helix)*cos(M_PI/2+theta_phi[i])*sin(theta_r[i]) + c_helix*cos(theta_r[i]);
        }

        if (!EVPA_rotation && !DopDep)
        {

            Proj_theta_B[0] = atan2(cos(t_helix + phi_helix),(c_helix*sin(deg2rad(theta_obs))/R - cos(deg2rad(theta_obs))*sin(t_helix + phi_helix))); //helical B in middle block
            for (i=1; i<(n_blocks); i++)
            {
                Proj_theta_B[i] = atan2(b_00[i],(c_00[i]*sin(deg2rad(theta_obs)) + a_00[i]*cos(deg2rad(theta_obs)))); //random B-fields in other blocks
            }

        }
        if (EVPA_rotation && !DopDep) //now starting a rotation, big chunk of blocks will be the helical field
        {
            for (l=0; l<(n_blocks/3-23); l++)
            { //1/3 is ratio of blocks which turn helical
                Proj_theta_B[l] = atan2(b_helical[l],(c_helical[l]*sin(deg2rad(theta_obs)) + cos(deg2rad(theta_obs))*a_helical[l])); //helical B in middle block
            }
            for (i=(n_blocks/3-23); i<(n_blocks); i++)
            {
                Proj_theta_B[i] = atan2(b_00[i],(c_00[i]*sin(deg2rad(theta_obs)) + a_00[i]*cos(deg2rad(theta_obs))));
            }

        } //Now each v vector is different for each block as theta_obs < theta_open
        if (!EVPA_rotation && DopDep)
        {
            Proj_theta_B[0] = DD_Beffective((c_helix*sin(deg2rad(theta_obs)) - R*cos(deg2rad(theta_obs))*sin(t_helix + phi_helix)), R*cos(t_helix + phi_helix), (R*sin(deg2rad(theta_obs))*sin(t_helix+phi_helix) + c_helix*cos(deg2rad(theta_obs))), cos(deg2rad(theta_obs))*sin(theta_r[i])*cos(theta_phi[i])+sin(deg2rad(theta_obs))*cos(theta_r[i]), sin(theta_r[i])*sin(theta_phi[i]), -sin(deg2rad(theta_obs))*sin(theta_r[i])*cos(theta_phi[i])+cos(deg2rad(theta_obs))*cos(theta_r[i]), gamma_bulk);
            for (i=1; i<(n_blocks); i++)
            {
                Proj_theta_B[i] = DD_Beffective((c_00[i]*sin(deg2rad(theta_obs)) + a_00[i]*cos(deg2rad(theta_obs))),b_00[i],(c_00[i]*cos(deg2rad(theta_obs))-a_00[i]*sin(deg2rad(theta_obs))),cos(deg2rad(theta_obs))*sin(theta_r[i])*cos(theta_phi[i])+sin(deg2rad(theta_obs))*cos(theta_r[i]), sin(theta_r[i])*sin(theta_phi[i]), -sin(deg2rad(theta_obs))*sin(theta_r[i])*cos(theta_phi[i])+cos(deg2rad(theta_obs))*cos(theta_r[i]),gamma_bulk); //random B-fields in other blocks
            }

        }
        if (EVPA_rotation && DopDep)
        {
            for (l=0; l<(n_blocks/3-23); l++)
            {
                Proj_theta_B[l] = DD_Beffective((c_helical[l]*sin(deg2rad(theta_obs)) + cos(deg2rad(theta_obs))*a_helical[l]), b_helical[l], (-sin(deg2rad(theta_obs))*a_helical[l] + c_helical[l]*cos(deg2rad(theta_obs))), cos(deg2rad(theta_obs))*sin(theta_r[i])*cos(theta_phi[i])+sin(deg2rad(theta_obs))*cos(theta_r[i]), sin(theta_r[i])*sin(theta_phi[i]), -sin(deg2rad(theta_obs))*sin(theta_r[i])*cos(theta_phi[i])+cos(deg2rad(theta_obs))*cos(theta_r[i]), gamma_bulk);
            }
            for (i=(n_blocks/3-23); i<(n_blocks); i++)
            {
                Proj_theta_B[i] = DD_Beffective((c_00[i]*sin(deg2rad(theta_obs)) + a_00[i]*cos(deg2rad(theta_obs))),b_00[i],(c_00[i]*cos(deg2rad(theta_obs))-a_00[i]*sin(deg2rad(theta_obs))),cos(deg2rad(theta_obs))*sin(theta_r[i])*cos(theta_phi[i])+sin(deg2rad(theta_obs))*cos(theta_r[i]), sin(theta_r[i])*sin(theta_phi[i]), -sin(deg2rad(theta_obs))*sin(theta_r[i])*cos(theta_phi[i])+cos(deg2rad(theta_obs))*cos(theta_r[i]),gamma_bulk); //random B-fields in other blocks
            }

        }

        //********************************************** ALP IC losses *******************************************//
        
        for (i=0; i<array_size; i++)
        {
            f_c_IC[i] = 0.0;
            
            Ne_etot += Ne_e[i];
            n_rad += urad_array[i] /(f_pol[i]*6.63E-34);
            f_c_IC[i] =gamma_e[i] * gamma_e[i] * f_pol[i];
        }
        
        for (i=0; i<array_size; i++) //1.79
        {
            Ps_per_m_IC[i] = 0.0;
            Ps_per_m_IC[i] = Ps_per_m[i] * 0.85 * (gamma_e[i] * gamma_e[i]) * (Ne_e[i] / Ne_etot) * (urad_array[i] /(f_pol[i]*6.63E-34))/n_rad; //* 1.79 //3.3
            IC_losses[i] = Ps_per_m_IC[i];
        }
        
        
        //**********************************************************************************************//
        
        
        
        //some code to estimate what dx should be
        dx_set = findminelement(Ne_e, array_size);//gets the lowest non_zero element. Higher energy electrons radiate more rapidly


        dx_R = 0.05*(R0+x*tan(deg2rad(theta_open_p)))/tan(deg2rad(theta_open_p)); //ensures Rnew <= 1.05 Rold
        dx_P = (Ne_e[dx_set]*(E_elecs[dx_set]-E_elecs[dx_set-1])*1.6E-19)/(Ps_per_m_elec[dx_set] + Ps_per_m_IC[dx_set]); //based on radiative losses
        //printf("dx_R, dx_P \t%.5e\t%.5e \n", dx_P, (Ne_e[dx_set]*(E_elecs[dx_set]-E_elecs[dx_set-1])*1.6E-19)/(Ps_per_m[dx_set]));
  
        if (dx_R < dx_P)
        {
            dx = dx_R;
        }

        if (dx_P < dx_R)
        {
            dx = dx_P;
        }

        //printf("dx_R dx_P dx %.5e\t%.5e\t%.5e \n", dx_R, dx_P, dx);

        
        
        //change values for next loop
        x += dx;
        R_prev = R;
        R = R_new(R0, deg2rad(theta_open_p), x);

        //correct synchrotron for full length + IC
        for (i=0; i<array_size; i++)
        {
            Ps_per_m[i] *= dx;
            Sync_losses[i] *= dx;
            Ps_per_m_IC[i] *= dx;             //ALP
            IC_losses[i] *= dx;              //ALP
            P_perp[i] *= dx;  //polarisation powers
            P_para[i] *= dx;
	        fprintf(jfile, "\t%.5e", j[i]);
            fprintf(freqfile, "\t%.5e", f_c[i]); //these are Doppler boosted later
            fprintf(opfile, "\t%.5e", k[i]);//*(R/R_prev)*(R/R_prev));
            fprintf(powfile, "\t%.5e", Ps_per_m[i]*pow(doppler_factor, 4.0));//*pow((doppler(beta_e[i], deg2rad(theta_obs))), 4.0));//*Ne_e[i]);
            fprintf(ICfile, "\t%.5e", Ps_per_m_IC[i]*pow(doppler_factor, 4.0)); //ALP
            fprintf(ICfreqfile, "\t%.5e",f_c_IC[i]); //ALP
            fprintf(Pperpfile, "\t%.5e", P_perp[i]);//*pow(doppler_factor, 4.0)); //polarisation power now doppler boosted in python file
            fprintf(Pparafile, "\t%.5e",P_para[i]);//*pow(doppler_factor, 4.0));
            //fprintf(Proj_Bfile,"\t%.5e",Proj_theta_B);
            //fprintf(betadata, "\t%.5e", beta_e[i]);

            if (i==array_size-1) //puts the outputs for the next jet section on a new line in files
            {
                fprintf(freqfile, "\n");
                fprintf(opfile, "\n");
        		fprintf(jfile, "\n");
                fprintf(powfile, "\n");
                fprintf(ICfile, "\n");
                fprintf(ICfreqfile, "\n");
                fprintf(Pperpfile, "\n");
                fprintf(Pparafile, "\n");
                //fprintf(Proj_Bfile, "\n");
            }
        }
        //same thing as above but now just to save Proj_Bs as this has a different size (nblocks)
        for(i=0; i<n_blocks; i++){
            //save Projected B field angle for each jet section and block here
            fprintf(Proj_Bfile,"\t%.5e",Proj_theta_B[i]);
            if (i==n_blocks-1){
                fprintf(Proj_Bfile,"\n");
            }
        }



        //compute the synchrotron losses + IC LOSSES
        for (i=0; i<array_size; i++)
        {
            if (i==0)
            {
                //Ne_losses[i]=0.0; //cannot drop to lower bins
                Ne_losses[i]=(Sync_losses[i] + IC_losses[i])*(3E8/3E8)/((E_elecs[i]-me_eV)*1.6E-19); //can't radiate anymore
            }
            else
            {
	      //Ne_losses[i] = (Sync_losses[i])*(3E8/3E8)/((E_elecs[i]-E_elecs[i-1])*1.6E-19);
	      Ne_losses[i] = (Sync_losses[i] + IC_losses[i])*(3E8/3E8)/((E_elecs[i]-E_elecs[i-1])*1.6E-19);
            }


            if (i==array_size-1)
            {
                Ne_gains[i] = 0.0;
            }
            else
            {
	      //Ne_gains[i] = (Sync_losses[i+1])*(3E8/3E8)/((E_elecs[i+1]-E_elecs[i])*1.6E-19);
	      Ne_gains[i] = (Sync_losses[i+1] + IC_losses[i+1])*(3E8/3E8)/((E_elecs[i+1]-E_elecs[i])*1.6E-19);
            }

	    u_rad = 0.0;
        }



        for (i=0; i<array_size; i++)
        {
            Ne_e[i] = Ne_e[i] - Ne_losses[i] + Ne_gains[i];


            if (x>4.9E12 && intermed_param==0) //arbitrary -> to replicate paper results
            {
                Ne_intermed[i] = Ne_e[i];

                if (i==array_size-1)
                {
                    intermed_param = 1; //ensures condition only met once
                }
            }

        }


        cleanpop(Ne_e, array_size, 10.0); //any elements with Ne_e<10 set to zero to avoid nans
        cleanpop(Ne_intermed, array_size, 10.0); //any elements with Ne_e<10 set to zero to avoid nans
        

        
	    B_prev = B;//need to update dfreqs
        B = get_newB(B0, R0, R);
        eps = epsilon(B);


        for (i=0; i<array_size; i++)
        {
            f_c[i] = f_crit(gamma_e[i], B);
            f_em_min[i] *= B/B_prev;//steps[nSteps];
            f_em_max[i] *= B/B_prev;//steps[nSteps];
            dfreqs[i] = f_em_max[i]-f_em_min[i];
            A_elecs[i] = A_orig*(Ne_e[i]/Ne_orig[i]);
            dN_dE[i] = electron_PL(A_elecs[i], alpha, E_elecs[i]*1.6E-19, E_max*1.6E-19);
            P_para[i] = 0.0;
            P_perp[i] = 0.0;


        }
        nSteps+=1; //allows the number of jet sections to be determined


        //save data needed at every jet section to solve line of sight opacity
        fprintf(basicdata, "\t%.5e\t%.5e\t%.5e\t%.5e \n", dx, x, B, R);

    }

    printf("The jet was divided into %d %s \n", nSteps, "section(s).");

    //store electron information
    for (i=0; i<array_size; i++)

    {
        fprintf(elecpop, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e \n", E_elecs[i], Ne_orig[i], Ne_intermed[i], Ne_e[i], dN_dE_orig[i]); //file elecpop.txt

    }

}


