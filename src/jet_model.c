/* 
* Program generating quiescent synchrotron and SSC (thomson limit) power spectrum and polarization
* for an optically thin multizone conical jet (B-field direction independent of x for SSC). 
* Outputs Total and individual Zone emission.
* File I/O built for parallelization on Slurm Cluster.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "jet_fns.h"
#include "mtwister.h" //random number generation
#include <time.h>
#include <omp.h>

#define Me_EV 0.511E6 //electron rest energy
#define C 3.0E8 //speed of light
#define Qe 1.6E-19 //electron charge
#define Me 9.11E-31 //electron mass kg
#define H 6.63E-34 //Plancks constant

const int ARRAY_SIZE = 50; // sets the number of synchrotron & SSC bins
const double ARRAY_SIZE_D = 50.0; // use to set log ratios, should be the same as above line but with .0
int i, l, m, n, o, p, g, h, nn; //some looping parameters
double c, d, q;
int dx_set; // use to define smallest non zero population

int main(int argc,char* argv[]) //argc is integer number of arguments passed, argv[0] is program name, argv[1..n] are arguments passed in string format
{
    // Jet Parameters
    double W_j = 1.3E37; // W jet power in lab frame. Should be OBSERVED POWER
    double L_jet = 5E20; // length in m in the fluid frame
    double E_min = 5.11E6; // Minimum electron energy 
    double E_max = 1.7E10; // Energy of the ECO in eV 
    double alpha = 1.85; // PL index of electrons
    double theta_open_p = 40.0; // opening angle of the jet in the fluid frame 
    double gamma_bulk = 16.0; // bulk Lorentz factor of jet material
    double B = 1E-4, B0 = 1E-4; // B-field at jet base
    double R0 = 0.0, R = 0.0;  // Radius of the jet at the base 3.32 works fairly well 
    double B_prev = 0.0; // changing parameters of the jet-initialise. R prev corrects for increasing jet volume 
    double theta_obs;  // observers angle to jet axis in rad 
    int N_BLOCKS; // for the TEMZ model, can have 1,7,19,37,61,91,127 blocks, (rings 0,1,2,3,4,5,6)
    int N_RINGS; // up to 6 rings possible atm, must choose number of rings corresponding to number of zones
    int SSC;
    sscanf(argv[4], "%lf", &theta_obs);
    sscanf(argv[5], "%d", &N_BLOCKS);
    sscanf(argv[6], "%d", &N_RINGS);
    sscanf(argv[10], "%d", &SSC);
    printf("Jet Parameters: %.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%d\n", W_j, gamma_bulk, theta_obs, theta_open_p, B0, E_max, alpha, N_BLOCKS);

    // Define radius at jet base
    R0 = R_0(W_j*(3.0/4.0), 1.0, gamma_bulk, B0);
    R = R_0(W_j*(3.0/4.0), 1.0, gamma_bulk, B0);
    double R_prev = R0; //initialize
    printf("Radius at jet base: %.5e \n", R);

    // Read in any required data (only Bessel Functions for now)
    FILE *FGfile;
    FGfile = fopen("src/FG.txt","r"); //Bessel Functions for Synchrotron

    double FG[75][3]; //Bessel Functions F, G and FG
    for (m=0; m<75; m++)
    {
           fscanf(FGfile, "%lf\t%lf\t%lf\n", &q, &c, &d);
           FG[m][0] = q;
           FG[m][1] = c;
           FG[m][2] = d;
    }

    // Output filenames for cluster
    char thread_idst[10];
    char task_idst[10];
    char results_dir[15];
    sscanf(argv[2], "%s", thread_idst);
    sscanf(argv[7], "%s", task_idst);
    sscanf(argv[9], "%s", results_dir);
    strcat(task_idst,"_");
    strcat(task_idst,thread_idst);
    strcat(task_idst,".txt");			
		
    char frange[40] = "/freqrange";
    strcat(frange,task_idst);
    prepend(frange, results_dir);
    char bdata[40] = "/basicdata";
    strcat(bdata,task_idst);
    prepend(bdata, results_dir);
    char kparams[40] = "/keyparams";
    strcat(kparams,task_idst);
    prepend(kparams, results_dir);
    //char Edensity[30] = "results/energydensity";
    //strcat(Edensity,task_idst);
    char testfil[40] = "/pi";
    strcat(testfil,task_idst);
    prepend(testfil, results_dir);
    char icz[40] = "/IC_Z";
    strcat(icz,task_idst);
    prepend(icz, results_dir);
    char sz[40] = "/S_Z";
    strcat(sz,task_idst);
    prepend(sz, results_dir);

		
    // Define some files to store output data
    FILE *freqrange, *basicdata, *keyparams, *PI, *IC_Z, *S_Z; //, *Xfile;

    IC_Z = fopen(icz, "w");
    S_Z = fopen(sz, "w");
    freqrange = fopen(frange, "w");//store frequency bin boundaries
    basicdata = fopen(bdata, "w"); //store B, x, R etc
    keyparams = fopen(kparams, "w"); //store Lj, gamma_bulk, theta_obs and theta_open
    PI = fopen(testfil,"w");
    //energydensity = fopen(Edensity, "w");
    //Proj_Bfile = fopen("Proj_Bfile.txt", "w"); //proj/ected B field onto plane of the sky for each section
    //block_thetafile = fopen("block_thetafile.txt","w"); //saves the angle to the line of sight of each block
    //Xfile = fopen("Xfile.txt","w");
    fprintf(keyparams, "\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%d\t%d \n", W_j, gamma_bulk, theta_obs, theta_open_p, alpha, B, E_max, N_BLOCKS, ARRAY_SIZE);

    //define some useful parameters to gauge progress
    int nSteps = 0; //how many sections the jet has broken down into
    double x = 0; //progress along the jet
    double dx = 0; //incremenet of step length
    double eps; //stores the eqn 2.18b
    double beta_bulk = lortovel(gamma_bulk);
    double doppler_factor = doppler(beta_bulk, deg2rad(theta_obs));
    int intermed_param = 0; //allows intermediate population to be determined to compare to paper
    clock_t time_end, time_begin; //to gauge timing for cluster
    double time_spent; //time spent on a single section calculation

    //define some parameters for determining the dx of each jet section
    double dx_R; //dx where R_new = 1.05*R_old
    double dx_P; //dx to ensure smallest bin population not depleted

    //some parameters used to define the population of electrons
    double A_elecs[ARRAY_SIZE]; //stores PL coefficient
    double Ne_e[ARRAY_SIZE]; // no of electrons in each bin
    double Ne_orig[ARRAY_SIZE]; // at jet base: allows comparison later
    double Ne_intermed[ARRAY_SIZE]; // intermediate population
    double dN_dE[ARRAY_SIZE]; // no of electrons per J in each bin
    double dN_dE_orig[ARRAY_SIZE];
    double dEe[ARRAY_SIZE]; //size of bin widths in eV
    double Ne_losses[ARRAY_SIZE]; //losses per section
    double Ne_gains[ARRAY_SIZE]; //gains per section
    double A_orig; // PL prefactor initially the same for all electrons
    double Ue_jetbase = 0.0; // energy in particles at conical region base
    double UB_jetbase = 0.0; //energy in particles at the jet base

    // parameters used to initialise electron population
    double Ee_max= 4.0E12;//Me_EV*3.0E4;//4.0E12; //sets the maximum electron energy in eV
    double ratio = pow((Ee_max/E_min), (1/ARRAY_SIZE_D)); //logarithmic ratio of electron bins: USED FOR BOUNDS
    double E_bounds[ARRAY_SIZE+1], indices[ARRAY_SIZE+1];
    arange(indices, ARRAY_SIZE+1); //as np.arange
    
    
    for (i=0; i<ARRAY_SIZE+1; i++)
    {
    E_bounds[i] = E_min*(pow(ratio, indices[i]));
    //printf("indices, bounds: %.5e\t%.5e \n", indices[i], E_bounds[i]); 
    }
    
    //now sort the bounds into max, min and E_elecs
    double E_elecs[ARRAY_SIZE]; // stores the energy in eV of each electron
    double E_elec_min[ARRAY_SIZE], E_elec_max[ARRAY_SIZE]; //these map onto critical frequencies, should be in eV
    
    for (i=0; i<ARRAY_SIZE; i++)
    {
    E_elec_min[i] = E_bounds[i];
    E_elec_max[i] = E_bounds[i+1];
    E_elecs[i] = pow(E_elec_min[i]*E_elec_max[i], 0.5); 
    }

    //parameters based on individual electron bins
    double gamma_e[ARRAY_SIZE]; // store Lorentz factors
    double beta_e[ARRAY_SIZE]; // v/c

    //parameters used in determining synchrotron emission
    double f_c[ARRAY_SIZE]; //store the critical frequencies
    double dfreqs[ARRAY_SIZE];// store the bin widths in critical frequency
    double f_em_min[ARRAY_SIZE], f_em_max[ARRAY_SIZE], f_em_min_IC[ARRAY_SIZE], f_em_max_IC[ARRAY_SIZE]; // store bin boundaries for critical frequencies
    double j[ARRAY_SIZE]; // emissivity per Hz per unit volume
    double k[ARRAY_SIZE]; // opacity, units of m^-1
    double Ps_per_m[ARRAY_SIZE]; // Synchrotron power per unit m
    double Ps_per_m_elec[ARRAY_SIZE]; //power lost by each electron bin per m
    double Ps_per_m_test[ARRAY_SIZE]; //used to normalise IC power emission
    memset(Ps_per_m_test, 0.0, ARRAY_SIZE*sizeof(Ps_per_m_test[0]));
    double Sync_losses[ARRAY_SIZE];// store total electron energy losses
    double P_single[ARRAY_SIZE];// power emitted by a single electron of energy E
    //double Ps_per_m_IC[ARRAY_SIZE]; //IC power per m
    double Ps_per_mIC_elec[ARRAY_SIZE]; //energy lost by each electron bin per m of IC emission
    memset(Ps_per_mIC_elec, 0.0, ARRAY_SIZE*sizeof(Ps_per_mIC_elec[0]));
    double IC_losses[ARRAY_SIZE]; //store electron IC energy losses

    //******************************** Inititalise seed photon energy density variables *************************************//

    //want to allocate these as dynamic arrays so can change their size with realloc
    double *R_array, *dx_array, **Pperp_array, **Ppara_array;
    double *ptri,**ptrii;
    double counter = 0.0;
    double dxsum = 0;
    double Area0;
    double dist = 0;
    
    //double urad_array_perpTEST[ARRAY_SIZE]; //photon energy density
    //memset(urad_array_perpTEST, 0.0, ARRAY_SIZE*sizeof(urad_array_perpTEST[0]));
    double urad_array_perp[N_BLOCKS][ARRAY_SIZE]; //photon energy density
    memset(urad_array_perp, 0.0, N_BLOCKS*ARRAY_SIZE*sizeof(urad_array_perp[0][0]));
    double urad_array_para[N_BLOCKS][ARRAY_SIZE];
    memset(urad_array_para, 0.0, N_BLOCKS*ARRAY_SIZE*sizeof(urad_array_para[0][0]));
    double dfactor_perp[N_BLOCKS][N_BLOCKS][ARRAY_SIZE]; //photon energy density
    memset(dfactor_perp, 0.0, N_BLOCKS*N_BLOCKS*ARRAY_SIZE*sizeof(dfactor_perp[0][0][0]));
    double dfactor_para[N_BLOCKS][N_BLOCKS][ARRAY_SIZE];
    memset(dfactor_para, 0.0, N_BLOCKS*N_BLOCKS*ARRAY_SIZE*sizeof(dfactor_para[0][0][0]));
    //double dfactor_perptest[N_BLOCKS][ARRAY_SIZE];
    //memset(dfactor_perptest, 0.0, N_BLOCKS*ARRAY_SIZE*sizeof(dfactor_perptest[0][0]));
    double dfactor_temp_perp = 0;
    double dfactor_temp_para = 0;
    double ang_factor; //reduction in power from small solid angle of zone already accounted for by dx/dxsum in dfactor__, this term remedies

    int buffer_subset[N_BLOCKS]; //different blocks have different buffer_sizes depending on their position in jet cross section
    for (i=0; i<N_BLOCKS; i++){ //cant allocate non zero with memset
        buffer_subset[i] = 1;
    }
    int buffer_size = 1;

    //****************Initialise Polarisation variables *****************************************************************************//
    //first set frequency bins for polarised powers to drop into:
    double f_pol[ARRAY_SIZE];
    double f_pol_IC[ARRAY_SIZE]; //fixed bins
    double f_pol_ICbounds[ARRAY_SIZE+1]; //for calculating IC bins properly
    double P_perp[ARRAY_SIZE];
    memset(P_perp, 0.0, ARRAY_SIZE*sizeof(P_perp[0]));
    double P_para[ARRAY_SIZE];
    memset(P_para, 0.0, ARRAY_SIZE*sizeof(P_para[0])); //initialise to make sure junk values arent inside upon first +=
    double P_perpIC[ARRAY_SIZE];
    memset(P_perpIC, 0.0, ARRAY_SIZE*sizeof(P_perpIC[0]));
    double P_paraIC[ARRAY_SIZE];
    memset(P_paraIC, 0.0, ARRAY_SIZE*sizeof(P_paraIC[0]));
    double P_X[ARRAY_SIZE];
    memset(P_X, 0.0, ARRAY_SIZE*sizeof(P_X[0]));
    double dfreqs_pol[ARRAY_SIZE]; //again fixed size
    double dfreqs_polIC[ARRAY_SIZE];
    double F_min;

    //cos(th_k) for compton integral
    double cosk[10];
    double phik[10];
    double dcosk = 1/5.5; //1/25.5;
    double dphik = (2*M_PI)/10;//(2*M_PI)/ARRAY_SIZE;
    for (i=0; i<10; i++){
        cosk[i] = (i-5)/5.5; //(i-25)/25.5;
        phik[i] = -M_PI + i*(2*M_PI)/10; //i*(2*M_PI)/ARRAY_SIZE;
    }


    //Fill electron arrays with initial parameters
    for (i=0; i<ARRAY_SIZE; i++)
    {
        E_elec_min[i] = E_elecs[i]*(1.0/pow(ratio, 0.5)); //lower electron eV bounds
        E_elec_max[i] = E_elecs[i]*pow(ratio, 0.5); //upper electron eV bounds
        dEe[i] = E_elec_max[i]-E_elec_min[i]; // sets bin widths in eV
        A_elecs[i] = A_PL(alpha, W_j, gamma_bulk, E_min*Qe, E_max*Qe, 1.0); //1 assumes equipartition. DIVIDED W_j by gamma^2 for jet frame
        dN_dE[i] = electron_PL(A_elecs[i], alpha, E_elecs[i]*Qe, E_max*Qe);
        dN_dE_orig[i] = electron_PL(A_elecs[i], alpha, E_elecs[i]*Qe, E_max*Qe);
        Ne_e[i] = electron_PL(A_elecs[i], alpha, E_elecs[i]*Qe, E_max*Qe)*dEe[i]*Qe;
        Ne_orig[i] = electron_PL(A_elecs[i], alpha, E_elecs[i]*Qe, E_max*Qe)*dEe[i]*Qe; //allows comparison of final to initial electron population

        //now calculate the initial critical frequencies
        gamma_e[i] = E_elec_min[i]/Me_EV; //works as all in eV
        f_em_min[i] = f_crit(gamma_e[i], B); //lower bin values for f_crits

        gamma_e[i] = E_elec_max[i]/Me_EV; //redefine for upper boundary
        f_em_max[i] = f_crit(gamma_e[i], B);

        gamma_e[i] = E_elecs[i]/Me_EV; //mid points
        beta_e[i] = lortovel(gamma_e[i]); // beta taken at midpts

        f_c[i] = f_crit(gamma_e[i], B); //assume all radiation occurs at the critical frequency
        dfreqs[i] = f_em_max[i]-f_em_min[i]; //obtain frequency bin widths

        f_pol[i] = f_c[i]; //fixed bin polarization equivalents of the above
        dfreqs_pol[i] = dfreqs[i];
        
        Ue_jetbase += Ne_e[i]*E_elecs[i]*Qe;
    }

    // SSC Frequency bins
    double bindiff;
    for (i=0; i<ARRAY_SIZE; i++) {
        bindiff = log10(gamma_e[1] * gamma_e[1] * (f_em_min[1]))-log10(gamma_e[0] * gamma_e[0] * (f_em_max[0]));
        f_em_min_IC[i] = pow(10,log10(gamma_e[i] * gamma_e[i] * f_em_min[i])-bindiff/2);
        f_em_max_IC[i] = pow(10,log10(gamma_e[i] * gamma_e[i] * f_em_max[i])+bindiff/2);
        dfreqs_polIC[i] = (f_em_max_IC[i]-f_em_min_IC[i]);
        f_pol_IC[i] = pow(10,log10(gamma_e[i] * gamma_e[i] * f_pol[i]));

        fprintf(freqrange,  "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e \n", f_em_min[i], f_em_max[i], f_c[i], dfreqs[i],
                f_em_min_IC[i], f_em_max_IC[i], f_pol_IC[i], dfreqs_polIC[i]);  //saves the frequency bins!
    }

    A_orig = A_elecs[0]; //all elements the same to begin with


    //print the energies to test equipartition
    UB_jetbase = R0*R0*C*B0*B0/(2*1E-7);
    //printf("Obs power, particles, B-fields: %.5e\t%.5e\t%.5e\t%.5e \n", W_j/(gamma_bulk*gamma_bulk), Ue_jetbase, UB_jetbase, Ue_jetbase+UB_jetbase);
    
    //redefine electron population to ensure equipartition
    for (i=0; i<ARRAY_SIZE; i++)
    {

        A_elecs[i] = A_PL(alpha, W_j, gamma_bulk, E_min*Qe, E_max*Qe, 1.0)*UB_jetbase/Ue_jetbase; //1 assumes equipartition. DIVIDED W_j by gamma^2 for jet frame
        dN_dE[i] = electron_PL(A_elecs[i], alpha, E_elecs[i]*Qe, E_max*Qe);
        Ne_e[i] = electron_PL(A_elecs[i], alpha, E_elecs[i]*Qe, E_max*Qe)*dEe[i]*Qe;
        Ne_orig[i] = electron_PL(A_elecs[i], alpha, E_elecs[i]*Qe, E_max*Qe)*dEe[i]*Qe; //allows comparison of final to initial electron population

    }
    
    //recalaculate equipartition fraction
    Ue_jetbase = 0.0; 
    for (i=0; i<ARRAY_SIZE; i++)
    {
        Ue_jetbase += Ne_e[i]*E_elecs[i]*Qe;
    }

    // set bins with fewer than ten electrons to zero: prevents resolution from being limited
    cleanpop(Ne_e, ARRAY_SIZE, 10.0);
    cleanpop(dN_dE, ARRAY_SIZE, 10.0);
    

    for (i=0; i<ARRAY_SIZE; i++)
    {
        A_elecs[i]*=(Ne_e[i]/Ne_orig[i]); //sets A to zero for tiny populations
    }
    
    //11/10/16 define some parameters to determine thick/thin jet sections
    double R_eff[ARRAY_SIZE]; //effective radius down to which can be seen -> smaller than jet radius if optically thick


    //HELICAL B-field: This variable will save the 3D B-fields on the sky for each jet zone, these B-fields evolve along the jet
    //becoming more transverse
    double Proj_theta_B[N_BLOCKS];
    //defining how far along the helix we are
    int helix_counter;
    sscanf(argv[2], "%d", &helix_counter); //taking input from command line how many 'days' we have been observing
    //parameters for helix
    double t_helix = helix_counter*M_PI/32; //helix parameter x=rcos(t+phi) y=rsin(t+phi) z =ct (r will be radius of jet)
    //Pi constant decided how quickly we sample along the helix -> timescale
    double phi_helix = 0;//M_PI/3; //phase shift, just an initial condition
    double c_helix = R0; //speed (ie number of coils per unit distance) -- high c => spaced out coils, low c => densely packed coils

    int EVPA_rotation;
    sscanf(argv[1], "%d", &EVPA_rotation); //taking input from command line if true begin EVPA rotation, if false do not
    int DopDep;
    sscanf(argv[3], "%d", &DopDep); //taking input from command line if true activate DopDep, if false do not

    //now the Bfield vectors are rotated about both the y and x axis each block with a different theta_obs
    //these give the angles each block is at with respect to the direction of the jet:
    //theta_x => rotation angle about x axis, theta_y => rotation angle about y axis
    double th = 2*atan(tan(deg2rad(theta_open_p))/gamma_bulk); //full jet diameter opening angle in lab frame in rad
    double theta_r[N_BLOCKS];                // Marscher configuration, rings of turbulent cells, middle 3 rings turbulent
    double theta_phi[N_BLOCKS];
    //TEMZ config, calculates theta_r and theta_phi for each block in the concentric circles
    for (i=0; i<N_BLOCKS; i++){
        if (i<1){
            theta_r[i] = 0;
            theta_phi[i] = 0;
        }
        else if (i>0 && i<7){

            theta_r[i] = th/(2*(N_RINGS+0.5));
            theta_phi[i] = -M_PI+(i-1)*2*M_PI/6;

        }
        else if (i>6 && i<19){

            theta_r[i] = 2*th/(2*(N_RINGS+0.5));
            theta_phi[i] = -M_PI+(i-7)*2*M_PI/12;

        }
        else if (i>18 && i<37){

            theta_r[i] = 3*th/(2*(N_RINGS+0.5));
            theta_phi[i] = -M_PI+(i-19)*2*M_PI/18;

        }
        else if (i>36 && i<61){

            theta_r[i] = 4*th/(2*(N_RINGS+0.5));
            theta_phi[i] = -M_PI+(i-37)*2*M_PI/24;

        }
        else if (i>60 && i<91){

            theta_r[i] = 5*th/(2*(N_RINGS+0.5));
            theta_phi[i] = -M_PI+(i-61)*2*M_PI/30;

        }
        else if (i>90 && i<127){

            theta_r[i] = 6*th/(2*(N_RINGS+0.5));
            theta_phi[i] = -M_PI+(i-91)*2*M_PI/36;

        }
    }

    double theta_tot[N_BLOCKS]; //total angle to line of sight of each block (for doppler boosting) in radians
    for (i=0; i<N_BLOCKS; i++){
        if (N_BLOCKS==1){
            theta_tot[i] = deg2rad(theta_obs);
        } else {
            theta_tot[i] = acos(cos(theta_r[i])*cos(deg2rad(theta_obs)) + sin(theta_r[i])*sin(deg2rad(theta_obs))*cos(M_PI - fabs(theta_phi[i])));
        }
    }

    //for (i=0; i<N_BLOCKS; i++){
    //    fprintf(block_thetafile, "\t%.5e", theta_tot[i]); //saving the angle to the direction of view of each block to be used in python file for doppler boost
    //}

    //alignment vector matrix for each block to every other block, for use in Synchrotron mixing between zones for IC emission, distances are in units of theta_r
    //assuming jet section is flat 2D (ok assumption given large section size) but have z component to rotate with theta_obs
    double align[N_BLOCKS][N_BLOCKS][3]; //vector from m zone to i zone jet frame
    double unitalign0[N_BLOCKS][N_BLOCKS][3]; //unit vector from m zone to i zone RPAR adjusted
    double unitalign1[N_BLOCKS][N_BLOCKS][3];// unit vector from m zone to i zone RPAR adjusted with offset cosk
    double unitalign2[N_BLOCKS][N_BLOCKS][3];
    double unitalign3[N_BLOCKS][N_BLOCKS][3];
    double unitalign4[N_BLOCKS][N_BLOCKS][3];
    double unitalign5[N_BLOCKS][N_BLOCKS][3];
    double unitalign6[N_BLOCKS][N_BLOCKS][3];
    for (i=0; i<N_BLOCKS; i++){
        for (m=0; m<N_BLOCKS; m++){
            align[i][m][0] = theta_r[i]*cos(theta_phi[i]) - theta_r[m]*cos(theta_phi[m]);//cos(deg2rad(theta_obs))*(theta_r[i]*cos(theta_phi[i]) - theta_r[m]*cos(theta_phi[m]));
            align[i][m][1] = theta_r[i]*sin(theta_phi[i]) - theta_r[m]*sin(theta_phi[m]);
            align[i][m][2] = 0.0;//-sin(deg2rad(theta_obs))*(theta_r[i]*cos(theta_phi[i]) - theta_r[m]*cos(theta_phi[m]));
        }
    }


    //Full and zonal Stokes vector bins for Synchrotron and IC respectively
    double S_Stokes[N_BLOCKS][ARRAY_SIZE][3];
    memset(S_Stokes, 0, sizeof(S_Stokes[0][0][0])* N_BLOCKS * ARRAY_SIZE * 3);
    double S_StokesTotal[ARRAY_SIZE][3];
    memset(S_StokesTotal, 0, sizeof(S_StokesTotal[0][0])* ARRAY_SIZE * 3);
    double S_Pi[ARRAY_SIZE];
    double S_PA[ARRAY_SIZE];
    double S_P[ARRAY_SIZE];
    double ICS_Pi[ARRAY_SIZE];
    double ICS_PA[ARRAY_SIZE];
    double ICS_P[ARRAY_SIZE]; //for the total S + IC, have to think about adding powers with different frequency bin sizes
    double ICS_StokesTotal[ARRAY_SIZE][3];
    memset(ICS_StokesTotal, 0, sizeof(ICS_StokesTotal[0][0])* ARRAY_SIZE * 3);

    double IC_Stokes[N_BLOCKS][ARRAY_SIZE][3];
    memset(IC_Stokes, 0, sizeof(IC_Stokes[0][0][0])* N_BLOCKS * ARRAY_SIZE * 3);
    double IC_StokesTotal[ARRAY_SIZE][3];
    memset(IC_StokesTotal, 0, sizeof(IC_StokesTotal[0][0])* ARRAY_SIZE * 3);

    double IC_StokesZTotal[N_BLOCKS][ARRAY_SIZE][3];
    memset(IC_StokesZTotal, 0, sizeof(IC_StokesZTotal[0][0][0])* N_BLOCKS * ARRAY_SIZE * 3);
    double S_StokesZTotal[N_BLOCKS][ARRAY_SIZE][3];
    memset(S_StokesZTotal, 0, sizeof(S_StokesZTotal[0][0][0])* N_BLOCKS * ARRAY_SIZE * 3);
    double IC_Pi[ARRAY_SIZE];
    double IC_PA[ARRAY_SIZE];
    double IC_P[ARRAY_SIZE];
    double zone_doppler;
    double q_theta;
    double zeta;
    double psi;
    int phik_start, phik_end, cosk_start, cosk_end;
    double cosk_single[7];
    double phik_single[7];
    int cosk_list[7];
    int phik_list[7];
    double u_rad = 0.0;
    double Blength; //for normalization
    double logdif;
    int binshiftIC;
    int binshiftS;
    double dfactor;

    double X = 0.0;
    double KN = 1.0; //klein nishina factor
    double v_k[3] = {0}; //incoming photon vector
    double e_k[3] = {0}; //incoming polarization vector perpendicular to B
    double _e[2] = {0}; //outgoing polarization vector
    double e_kpara[3] = {0}; //incoming polarization vector parallel to B
    double zone_d; //distance from other zones
    //now including full RPAR treatment for IC
    double Btheta = 90*M_PI/180;
    double Pperpperp = 0.0;
    double Pperppara = 0.0;
    double Pparaperp = 0.0;
    double Pparapara = 0.0;
    double effective_alpha[ARRAY_SIZE];
    
    int breakstep;
    sscanf(argv[8], "%d", &breakstep);

    //choosing random vectors (B_0,B_1,B_2) in unit sphere for random blocks initial B directions:
    int thread_id;
    int task_id;
    sscanf(argv[2], "%d", &thread_id);
    sscanf(argv[7], "%d", &task_id);
    MTRand seedr = seedRand((unsigned)time(NULL)+(unsigned)(1631*task_id + 335*thread_id)); //random seed supplemented by task and thread id
    //MTRand seedr = seedRand(11); //fix random seed
    int nLT_sections = 30;

    double Bx_helical[N_BLOCKS]; //simplifies angular shifts for the helical components
    double By_helical[N_BLOCKS];
    double Bz_helical[N_BLOCKS];
    double BX[nLT_sections][N_BLOCKS]; //for joint helical and random array
    double BY[nLT_sections][N_BLOCKS];
    double BZ[nLT_sections][N_BLOCKS];
    //B_effectives due to RPAR for each zone in every other zone
    double B_effectives[nLT_sections][N_BLOCKS][N_BLOCKS][3];
    memset(B_effectives, 0, sizeof(B_effectives[0][0][0][0])* nLT_sections *  N_BLOCKS * N_BLOCKS * 3);

    int marker = 0;
    int marker_prev = 0;
    int marker_list[N_BLOCKS];
    for (n=0; n<N_BLOCKS; n++) {
        marker_list[n] = nLT_sections / 2;
        printf("marker_list init = %d", marker_list[n]);
    }

    double B_0[nLT_sections][N_BLOCKS];
    double B_1[nLT_sections][N_BLOCKS];
    double B_2[nLT_sections][N_BLOCKS];
    double unif_theta;
    double unif_phi;
    for (l=0; l<nLT_sections; l++){
        for (i=0; i<(N_BLOCKS); i++){ //these are the random B-field vectors in each block as if we are looking straight down jet (have to rotate by theta_obs)
          unif_phi = 2 * M_PI * genRand(&seedr);
          unif_theta = acos(2 * genRand(&seedr) - 1);

          B_0[l][i] = cos(unif_phi) * sin(unif_theta);
          B_1[l][i] = sin(unif_phi) * sin(unif_theta);
          B_2[l][i] = cos(unif_theta);

          //printf("randBs %.5e\t%.5e\t%.5e\n",B_0[i],B_1[i],B_2[i]);
       }
    }

    //Putting B_effectives setup outside of emission loop because the B-fields do not change along the jet length for ICPaper
    //Rotate B-fields along spherical cone surface slightly and also transversify depending on R/R0 (not for IC paper)
    for (i=0; i<(N_BLOCKS); i++){
          Bx_helical[i] = -sin(t_helix+phi_helix) / sqrt(2);
          By_helical[i] = cos(t_helix+phi_helix) / sqrt(2); //no transversification or cone surface, 45 deg helix
          Bz_helical[i] = 1 / sqrt(2);
    }

    if (!EVPA_rotation) {
        for (l=0; l<nLT_sections; l++){
             for (i=0; i<(N_BLOCKS); i++) {
                     BX[l][i] = (B_2[l][i]*sin(deg2rad(theta_obs)) + B_0[l][i]*cos(deg2rad(theta_obs)));
                     BY[l][i] = B_1[l][i];
                     BZ[l][i] = (B_2[l][i]*cos(deg2rad(theta_obs))-B_0[l][i]*sin(deg2rad(theta_obs)));
                     printf("BXBYBZ %d\t%d\t%.5e\t%.5e\t%.5e\n",l,i,BX[l][i],BY[l][i],BZ[l][i]);
             }
        }
    } else {
        for (i=0; i<(N_BLOCKS); i++) {
            if (i<19){
                BX[l][i] = (Bz_helical[i]*sin(deg2rad(theta_obs)) + Bx_helical[i]*cos(deg2rad(theta_obs)));
                BY[l][i] = By_helical[i];
                BZ[l][i] = (Bz_helical[i]*cos(deg2rad(theta_obs))-Bx_helical[i]*sin(deg2rad(theta_obs)));
            } else {
                BX[l][i] = (B_2[l][i]*sin(deg2rad(theta_obs)) + B_0[l][i]*cos(deg2rad(theta_obs)));
                BY[l][i] = B_1[l][i];
                BZ[l][i] = (B_2[l][i]*cos(deg2rad(theta_obs))-B_0[l][i]*sin(deg2rad(theta_obs)));
            }

        }
    }

    if (!DopDep){
        for (i=0; i<(N_BLOCKS); i++) {
            for (l=0; l<N_BLOCKS; l++){
                Blength = sqrt(BX[m][i]*BX[m][i]+BY[m][i]*BY[m][i]+BZ[m][i]*BZ[m][i]); //normalization important to get correct z angle

                B_effectives[m][l][i][0] = BX[m][i]/Blength;
                B_effectives[m][l][i][1] = BY[m][i]/Blength;
                B_effectives[m][l][i][2] = BZ[m][i]/Blength;

            }
        }
    } else {
        for (i=0; i<N_BLOCKS; i++){
            for (l=0; l<(N_BLOCKS); l++) {

                for (m=0; m<nLT_sections; m++){
                    DD_3Beffective(BX[m][i],BY[m][i],BZ[m][i],
                                cos(deg2rad(theta_obs))*sin(theta_r[l])*cos(theta_phi[l])+sin(deg2rad(theta_obs))*cos(theta_r[l]),
                                sin(theta_r[l])*sin(theta_phi[l]), -sin(deg2rad(theta_obs))*sin(theta_r[l])*cos(theta_phi[l])+cos(deg2rad(theta_obs))*cos(theta_r[l]),
                                gamma_bulk,B_effectives[m][l][i]);
                    }

                 DD_3Beffective(align[l][i][0],align[l][i][1],align[l][i][2],
                             cos(deg2rad(theta_obs))*sin(theta_r[l])*cos(theta_phi[l])+sin(deg2rad(theta_obs))*cos(theta_r[l]),
                             sin(theta_r[l])*sin(theta_phi[l]), -sin(deg2rad(theta_obs))*sin(theta_r[l])*cos(theta_phi[l])+cos(deg2rad(theta_obs))*cos(theta_r[l]),
                             gamma_bulk,unitalign0[l][i]);

                 DD_3Beffective(align[l][i][0],align[l][i][1], sin(18.2*M_PI/180) * sqrt(pow(align[l][i][0],2) + pow(align[l][i][1],2)), //adding cosk offsets
                             cos(deg2rad(theta_obs))*sin(theta_r[l])*cos(theta_phi[l])+sin(deg2rad(theta_obs))*cos(theta_r[l]),
                             sin(theta_r[l])*sin(theta_phi[l]), -sin(deg2rad(theta_obs))*sin(theta_r[l])*cos(theta_phi[l])+cos(deg2rad(theta_obs))*cos(theta_r[l]),
                             gamma_bulk,unitalign1[l][i]);
                 
                 DD_3Beffective(align[l][i][0],align[l][i][1], sin(36.4*M_PI/180) * sqrt(pow(align[l][i][0],2) + pow(align[l][i][1],2)),
                             cos(deg2rad(theta_obs))*sin(theta_r[l])*cos(theta_phi[l])+sin(deg2rad(theta_obs))*cos(theta_r[l]),
                             sin(theta_r[l])*sin(theta_phi[l]), -sin(deg2rad(theta_obs))*sin(theta_r[l])*cos(theta_phi[l])+cos(deg2rad(theta_obs))*cos(theta_r[l]),
                             gamma_bulk,unitalign2[l][i]);

                 DD_3Beffective(align[l][i][0],align[l][i][1], sin(-18.2*M_PI/180) * sqrt(pow(align[l][i][0],2) + pow(align[l][i][1],2)),
                             cos(deg2rad(theta_obs))*sin(theta_r[l])*cos(theta_phi[l])+sin(deg2rad(theta_obs))*cos(theta_r[l]),
                             sin(theta_r[l])*sin(theta_phi[l]), -sin(deg2rad(theta_obs))*sin(theta_r[l])*cos(theta_phi[l])+cos(deg2rad(theta_obs))*cos(theta_r[l]),
                             gamma_bulk,unitalign3[l][i]);
                 
                 DD_3Beffective(align[l][i][0],align[l][i][1], sin(-36.4*M_PI/180) * sqrt(pow(align[l][i][0],2) + pow(align[l][i][1],2)),
                             cos(deg2rad(theta_obs))*sin(theta_r[l])*cos(theta_phi[l])+sin(deg2rad(theta_obs))*cos(theta_r[l]),
                             sin(theta_r[l])*sin(theta_phi[l]), -sin(deg2rad(theta_obs))*sin(theta_r[l])*cos(theta_phi[l])+cos(deg2rad(theta_obs))*cos(theta_r[l]),
                             gamma_bulk,unitalign4[l][i]);

                 DD_3Beffective(align[l][i][0],align[l][i][1], sin(-60*M_PI/180) * sqrt(pow(align[l][i][0],2) + pow(align[l][i][1],2)),
                             cos(deg2rad(theta_obs))*sin(theta_r[l])*cos(theta_phi[l])+sin(deg2rad(theta_obs))*cos(theta_r[l]),
                             sin(theta_r[l])*sin(theta_phi[l]), -sin(deg2rad(theta_obs))*sin(theta_r[l])*cos(theta_phi[l])+cos(deg2rad(theta_obs))*cos(theta_r[l]),
                             gamma_bulk,unitalign5[l][i]);
                 
                 DD_3Beffective(align[l][i][0],align[l][i][1], sin(60*M_PI/180) * sqrt(pow(align[l][i][0],2) + pow(align[l][i][1],2)),
                             cos(deg2rad(theta_obs))*sin(theta_r[l])*cos(theta_phi[l])+sin(deg2rad(theta_obs))*cos(theta_r[l]),
                             sin(theta_r[l])*sin(theta_phi[l]), -sin(deg2rad(theta_obs))*sin(theta_r[l])*cos(theta_phi[l])+cos(deg2rad(theta_obs))*cos(theta_r[l]),
                             gamma_bulk,unitalign6[l][i]);
                
            }
        }
    }

//        for (i=0; i<(N_BLOCKS); i++) { //save B-fields for visualization
//            fprintf(Proj_Bfile,"\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n",theta_r[i],theta_phi[i],
//                    theta_tot[i],atan2(B_1[i],B_0[i]),atan2(By[i],Bx[i]),atan2(By_helical[i],Bx_helical[i]),
//                    atan2(BY[i],BX[i]),atan2(B_effectives[i][i][1],B_effectives[i][i][0]));
//        }


    // Actual Jet emission calculation loop
    printf("Beginning jet analysis... \n");
    while (x<L_jet)
    {
        time_begin = clock();
        eps = epsilon(B);//same for all populations

	    //----------------- B-field projection onto sky for this section --------------//
        // Now have introduced R dependence on the random blocks as well, works just like for the helical case, just * R_old/R_new on B_2 term
        //now can also include doppler depolarisation, changing the B field to EFFECTIVE B fields which give correct perp E component as if it had
        //been doppler depolarisedx
        //Each block has its own theta_obs as we are looking inside the conical jet

        //Jet section mixing here: if Rtan(theta_obs) > R_0/sqrt(N) then we start assigning new B-fields to edge zones and so on, equivalent to seeing emission from zones in section
        //behind and in front.
        if (x > 2 * R0 / (sqrt(N_BLOCKS)) / (1 - beta_bulk)){ //test to make sure mixing due to finite gamma doesnt occur
            printf("\n WARNING: GAMMA MIXING! \n");
        }

        //******************************************************************************************************//
        //Light travel angular effect + gamma bulk effect: causes section mixing
        if (R*tan(deg2rad(theta_obs)) > R0 / (sqrt(N_BLOCKS)) ){
            for (i=0; i<(N_BLOCKS); i++){
                if (fabs(R * 2 * theta_r[i]*cos(theta_phi[i]) / th) > (abs(2*(marker_list[i]-14) - 1) * R0 / (sqrt(N_BLOCKS)) ) / tan(deg2rad(theta_obs))  ) {
                    if (fabs(theta_phi[i]) > M_PI/2) {
                       marker_list[i] += 1;
                    } else {
                       if (marker_list[i] == 15) {
                          marker_list[i] -= 2;
                       } else {
                         marker_list[i] -= 1;
                       }
                    }
                }
                //printf("i, zone_marker = %d \t  %d \n", i, marker_list[i]);
            }
           printf("\n WARNING: Section Mixing has begun! \n");
        }
        /***************************************************************************************************/

        //begin with Synchrotron emission
        for (i=0; i<ARRAY_SIZE; i++)
        {
            P_single[i] = larmor(B, beta_e[i], gamma_e[i]); //gives the radiated power by one electron in each bin
            j[i] = j_per_hz(P_single[i], R, dN_dE[i], eps, f_c[i]);//*(R_next/R)*(R_next/R);//gives emissivity per unit frequency

            k[i] = k_new(j[i], eps, f_c[i]); //gives the corresponding opacity
            Ps_per_m_elec[i] = P_single[i]*dN_dE[i]*E_elecs[i]*Qe*dfreqs_pol[i] / (2*f_c[i]*C);
            //determine whether jet is optically thin or thick -> this part is incomplete
            R_eff[i] = power_emitted(j[i], k[i], R)/j[i]; //depth down to which can be seen in one second

            //Polarisation -> Power emitted parallel and perpendicular to projected B field at each frequency interval per m (/ l_c) 
            for (l=0; l<ARRAY_SIZE; l++) {
                for (m=0; m<75; m++) {
                    if ((f_pol[l]/f_c[i]) >= FG[m][0] && (f_pol[l]/f_c[i]) <= FG[m+1][0]) {
                        P_perp[l] += (sqrt(3)*M_PI/(16*Me*pow(C,3)))*pow(Qe,3) * B * FG[m+1][1] * dN_dE[i]*dEe[i]*Qe * /*dfreqs_pol[l] */ 4.4216E11 / (1E-7);  //power per unit frequency, this is Power(freq), multiply by freq later to get nuF(nu)
                        P_para[l] += (sqrt(3)*M_PI/(16*Me*pow(C,3)))*pow(Qe,3) * B * FG[m+1][2] * dN_dE[i]*dEe[i]*Qe * /*dfreqs_pol[l] */ 4.4216E11 / (1E-7); //*dx !!! (this happens further down) and * other constants, already have /l_c
                    }

                    else if ((f_pol[l]/f_c[i]) <= FG[0][0] || (f_pol[l]/f_c[i]) >= FG[73][0]) {
                        P_perp[l] += 0.0;
                        P_para[l] += 0.0;
                    }
                }
            }
            Sync_losses[i] = Ps_per_m_elec[i];
            //printf("Ps_per_me_test %.5e\t%.5e\n", Ps_per_m_test[i],Ps_per_m_elec[i]);
        }


        // Find marker for LT section mixing
        marker = (int)round(R*tan(deg2rad(theta_obs)) /  (2 * R0 / (sqrt(N_BLOCKS)) ) );
	printf("MARKER: %d \n", marker);

        //some code to estimate what dx early for energy density calcultion, based off sync losses only, ok if IC losses not too big, early jet
        dx_set = findminelementNO0(Ne_e, ARRAY_SIZE);//gets the lowest non_zero element. Higher energy electrons radiate more rapidly

        dx_R = 0.05*(R0+x*tan(deg2rad(theta_open_p)))/tan(deg2rad(theta_open_p)); //ensures Rnew <= 1.05 Rold
        dx_P = (Ne_e[dx_set]*(E_elecs[dx_set]-E_elecs[dx_set-1])*Qe)/(Ps_per_m_elec[dx_set]); //based on radiative losses

        if (dx_R < dx_P){
           dx = dx_R;
        } else {
          dx = dx_P;
        }

	/*************************************************************************************************************/
        // Buffer of previous synchrotron emission setup
        if (x == 0) {
            // max buffer_size is when x == 2R_array[min]
            R_array = (double *)malloc(sizeof(double));
            dx_array = (double *)malloc(sizeof(double));
            Pperp_array = (double **)malloc(1 * sizeof(double *));
            Pperp_array[0] = (double *)malloc(ARRAY_SIZE * sizeof(double)); //dynamically allocating 2D double arrays
            Ppara_array = (double **)malloc(1 * sizeof(double *));
            Ppara_array[0] = (double *)malloc(ARRAY_SIZE * sizeof(double));

            if (R_array == NULL || dx_array == NULL) {
                printf("malloc failed\n");
                fprintf(stderr, "malloc failed\n");
                return(-1);
            }

            R_array[0] = R0;
            dx_array[0] = dx;
            for (l=0; l<ARRAY_SIZE; l++){
                Ppara_array[0][l] = P_para[l];
                Pperp_array[0][l] = P_perp[l];
            }

        } else if (2*theta_r[N_BLOCKS - 1]/th * R_array[0] + R_array[findminelement(R_array,buffer_size)] > x + dx) { //add rows to matrices, else do nothing
            for (n=0; n<N_BLOCKS; n++){
                if (2*theta_r[n]/th * R_array[0] + R_array[findminelement(R_array,buffer_size)] > x + dx){
                    buffer_subset[n] += 1;
                }
            }
            buffer_size += 1;

            ptri = (double *)realloc(R_array, buffer_size * sizeof(R_array)); // add one row
            R_array = ptri;

            ptri = (double *)realloc(dx_array, buffer_size * sizeof(dx_array)); // add one row
            dx_array = ptri;

            ptrii = (double **)realloc(Ppara_array, buffer_size * sizeof(*Ppara_array)); //increment rows by 1
            Ppara_array = ptrii;
            Ppara_array[buffer_size - 1] = (double *)malloc(ARRAY_SIZE * sizeof(double)); //allocate memory for new row

            ptrii = (double **)realloc(Pperp_array, buffer_size * sizeof(*Pperp_array)); //increment rows by 1
            Pperp_array = ptrii;
            Pperp_array[buffer_size - 1] = (double *)malloc(ARRAY_SIZE * sizeof(double)); //allocate memory for new row

        }

        //now shift rows in buffer down by 1, and assign new buffer values to top row afterwards
        if (x != 0 && buffer_size > 1) {
            for (n = buffer_size - 2; n>=0; n--){
                R_array[n+1] = R_array[n];
                dx_array[n+1] = dx_array[n];
                for (l=0; l<ARRAY_SIZE; l++){
                    Pperp_array[n+1][l] = Pperp_array[n][l];
                    Ppara_array[n+1][l] = Ppara_array[n][l];
                }
            }

            R_array[0] = R;
            dx_array[0] = dx;
            for (l=0; l<ARRAY_SIZE; l++){
                    Pperp_array[0][l] = P_perp[l];
                    Ppara_array[0][l] = P_para[l];
            }
        }

	/******************************************************************************************************************/
        Area0 = circle_area(2*theta_r[0]/th,dx_array[0],dx_array[0],R_array[0]);

        // Calculate effective alpha for electron population as function of energy and Calculate Synchrotron Stokes parameters 
        for (l=0; l<ARRAY_SIZE; l++){
            //P_para[l] = 0.2*P_perp[l];
            //make alpha dependent on frequency here
            if (l == 0){
               effective_alpha[l] = -(log10(dN_dE[l+1]) - log10(dN_dE[l]))/(log10(E_elecs[l+1]*Qe) - log10(E_elecs[l]*Qe));
            } else if (l > ARRAY_SIZE - 2) {
                effective_alpha[l] = -(log10(dN_dE[l]) - log10(dN_dE[l-1]))/(log10(E_elecs[l]*Qe) - log10(E_elecs[l-1]*Qe));
            } else {
                effective_alpha[l] = -(log10(dN_dE[l+1]) - log10(dN_dE[l-1]))/(log10(E_elecs[l+1]*Qe) - log10(E_elecs[l-1]*Qe));
            }
            if (isnan(effective_alpha[l]) || isinf(effective_alpha[l])) {
                effective_alpha[l] = 8.41283e+01;
            }
            //printf("alphaef \t%.5e\n",effective_alpha[l]);
            for (g=0; g<N_BLOCKS; g++) { //Sync Stokes Parameters
                q_theta = pow(sin(acos(B_effectives[marker_list[g]][g][g][2])),(effective_alpha[l]+1)/2); //factor to account for weaker emission when B_field pointed closer to line of sight
                S_Stokes[g][l][0] += q_theta*P_perp[l]/N_BLOCKS, S_Stokes[g][l][1] += (q_theta*P_perp[l]/N_BLOCKS)*cos(2*atan2(B_effectives[marker_list[g]][g][g][0],-B_effectives[marker_list[g]][g][g][1])), S_Stokes[g][l][2] += (q_theta*P_perp[l]/N_BLOCKS)*sin(2*atan2(B_effectives[marker_list[g]][g][g][0],-B_effectives[marker_list[g]][g][g][1]));
                S_Stokes[g][l][0] += q_theta*P_para[l]/N_BLOCKS, S_Stokes[g][l][1] += (q_theta*P_para[l]/N_BLOCKS)*cos(2*atan2(B_effectives[marker_list[g]][g][g][1],B_effectives[marker_list[g]][g][g][0])), S_Stokes[g][l][2] += (q_theta*P_para[l]/N_BLOCKS)*sin(2*atan2(B_effectives[marker_list[g]][g][g][1],B_effectives[marker_list[g]][g][g][0]));
                //printf("SSStokes123 %.5e\t%.5e\t%.5e\n", S_Stokes[g][l][0],S_Stokes[g][l][1],S_Stokes[g][l][2]);
            }


            //Energy Density loop, calculates the total energy density seen at each zone for this section of the jet
            // Always need to have uradARRAY[] = sum (dfactor over second index)
            // dfactor array should be all that is needed, dfactor[g][h][l] takes over role of urad_array_perp[l]
            // dfactor[g][h][l] = synchrotron energy density at zone h from zone g at frequency l
            for (n=0; n<N_BLOCKS; n++) {
                dxsum = 0;
                for (i=0; i<buffer_subset[n]; i++){
		            if (dxsum + dx_array[i] < 2*theta_r[n]/th * R_array[0] + R_array[i]){ //this takes place of buffer_subset, making sure emission from further than R_i not contributing
			            dxsum += dx_array[i];
                        urad_array_perp[n][l] += (Pperp_array[i][l] ) * dfreqs_pol[l] * dx_array[i] / (M_PI*pow(R_array[i],2)*C) * (circle_area(2*theta_r[n]/th * R_array[0],dx_array[i],dxsum, R_array[i]) / Area0) * (dx_array[0] / dxsum);
                        urad_array_para[n][l] += (Ppara_array[i][l] ) * dfreqs_pol[l] * dx_array[i] / (M_PI*pow(R_array[i],2)*C) * (circle_area(2*theta_r[n]/th * R_array[0],dx_array[i],dxsum, R_array[i]) / Area0) * (dx_array[0] / dxsum);
		    			
                      counter = 0.0; //divides up number zones contributing at the same time
                      for (g=0; g<N_BLOCKS; g++) {
                          dist = fabs( sqrt( pow( 2*align[g][0][0]/th*R_array[i] - 2*theta_r[n]/th*R_array[0]*cos(theta_phi[n]),2) + pow( 2*align[g][0][1]/th*R_array[i] - 2*theta_r[n]/th*R_array[0]*sin(theta_phi[n]),2) ) - dxsum );
                          if (dist < 1/(2*(N_RINGS + 0.5))*R_array[i]) {
                             counter += 1;
                          }
                      }

                      for (g=0; g<N_BLOCKS; g++) {
                          dist = fabs( sqrt( pow( 2*align[g][0][0]/th*R_array[i] - 2*theta_r[n]/th*R_array[0]*cos(theta_phi[n]),2) + pow( 2*align[g][0][1]/th*R_array[i] - 2*theta_r[n]/th*R_array[0]*sin(theta_phi[n]),2) ) - dxsum );
                          if (dist < 1/(2*(N_RINGS + 0.5))*R_array[i]){
                              dfactor_perp[n][g][l] += (Pperp_array[i][l] ) * dfreqs_pol[l] * dx_array[i] / (M_PI*pow(R_array[i],2)*C) * (circle_area(2*theta_r[n]/th * R_array[0],dx_array[i],dxsum, R_array[i]) / Area0) * (dx_array[0] / dxsum) / counter;
                              dfactor_para[n][g][l] += (Ppara_array[i][l] ) * dfreqs_pol[l] * dx_array[i] / (M_PI*pow(R_array[i],2)*C) * (circle_area(2*theta_r[n]/th * R_array[0],dx_array[i],dxsum, R_array[i]) / Area0) * (dx_array[0] / dxsum) / counter;
                          }
                      }
                    }
                }
            }
        }

        //fprintf(energydensity, "%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n", urad_array_perp[0][5],urad_array_perp[0][20],urad_array_perp[1][40],urad_array_perp[1][5],urad_array_perp[1][20],urad_array_perp[1][40],dfactor_perp[0][0][5]/urad_array_perp[0][5], dfactor_perp[1][1][5]/urad_array_perp[1][5], urad_array_perpTEST[5], urad_array_perpTEST[20], urad_array_perpTEST[40]);

        for (l=0; l<ARRAY_SIZE; l++) { //reset these values to be calculated this step
            P_perpIC[l] = 0.0;
            P_paraIC[l] = 0.0;
        }


        //********************************************** ALP SSC losses *******************************************//
        //XXX: given energy density prescription do we need 1/nblocks in Pperpperp etc? Yes, still do because Pperp_array uses full electron pop power

        if (SSC){

            #pragma omp parallel for private(h,n,l,p,m,i,o) reduction(+:IC_Stokes[:N_BLOCKS][:ARRAY_SIZE][:3])
            for (g=0; g<N_BLOCKS; g++){ //B-fields
                printf("NUM_THREADS: %d\n", omp_get_num_threads());
                double Btheta = 90*M_PI/180;
                double q_theta;
                double zeta;
                double cosk_single[7];
                double phik_single[7];
                int phik_start, phik_end, cosk_start, cosk_end;
                int cosk_list[7];
                int phik_list[7];
                int nn;
                double F_min;
                double KN = 1.0; //empirical klein nishina factor
                double v_k[3] = {0}; //incoming photon vector
                double e_k[3] = {0}; //incoming polarization vector perpendicular to B
                double _e[2] = {0}; //outgoing polarization vector
                double e_kpara[3] = {0}; //incoming polarization vector parallel to B
                double zone_d; //distance from other zones
                double Pperpperp = 0.0;
                double Pperppara = 0.0;
                double Pparaperp = 0.0;
                double Pparapara = 0.0;

                for (h=0; h<N_BLOCKS; h++){ //block locations
                    Btheta = acos(B_effectives[marker_list[h]][h][g][2]); //compton polarization fraction very dependent on this for a single zone, 90deg gives highest
                    zeta = atan(B_effectives[marker_list[h]][h][g][1]/B_effectives[marker_list[h]][h][g][0]); //to rotate each Stokes to lab frame


                    if (h==g){
                        phik_start = 0, phik_end = 10, cosk_start = 0, cosk_end = 10;
                        //rotate to B_effective
                        ang_factor = 1; //adjust so total IC power always the same

                    } else if (h!=g && dfactor_perp[h][g][0] != 0) {
                        //first rotate align vector in exactly the same rotation as B -> Beffective, taken care of by unitalign
                        //then rotate unitalign vector by -zeta about z axis to get B in x-z plane
                        cosk_single[0] = unitalign0[h][g][2]; //z component isnt affected by -zeta rotation
                        phik_single[0] = atan2(unitalign0[h][g][0]*sin(-zeta) + unitalign0[h][g][1]*cos(-zeta), unitalign0[h][g][0]*cos(-zeta) - unitalign0[h][g][1]*sin(-zeta)); //this is between -pi and pi, but phik between 0 and 2pi
                        cosk_single[1] = unitalign1[h][g][2]; //z component isnt affected by -zeta rotation
                        phik_single[1] = atan2(unitalign1[h][g][0]*sin(-zeta) + unitalign1[h][g][1]*cos(-zeta), unitalign1[h][g][0]*cos(-zeta) - unitalign1[h][g][1]*sin(-zeta));
                        cosk_single[2] = unitalign2[h][g][2]; //z component isnt affected by -zeta rotation
                        phik_single[2] = atan2(unitalign2[h][g][0]*sin(-zeta) + unitalign2[h][g][1]*cos(-zeta), unitalign2[h][g][0]*cos(-zeta) - unitalign2[h][g][1]*sin(-zeta));
                        cosk_single[3] = unitalign3[h][g][2]; //z component isnt affected by -zeta rotation
                        phik_single[3] = atan2(unitalign3[h][g][0]*sin(-zeta) + unitalign3[h][g][1]*cos(-zeta), unitalign3[h][g][0]*cos(-zeta) - unitalign3[h][g][1]*sin(-zeta));
                        cosk_single[4] = unitalign4[h][g][2]; //z component isnt affected by -zeta rotation
                        phik_single[4] = atan2(unitalign4[h][g][0]*sin(-zeta) + unitalign4[h][g][1]*cos(-zeta), unitalign4[h][g][0]*cos(-zeta) - unitalign4[h][g][1]*sin(-zeta));
                        cosk_single[5] = unitalign5[h][g][2]; //z component isnt affected by -zeta rotation
                        phik_single[5] = atan2(unitalign5[h][g][0]*sin(-zeta) + unitalign5[h][g][1]*cos(-zeta), unitalign5[h][g][0]*cos(-zeta) - unitalign5[h][g][1]*sin(-zeta));
                        cosk_single[6] = unitalign6[h][g][2]; //z component isnt affected by -zeta rotation
                        phik_single[6] = atan2(unitalign6[h][g][0]*sin(-zeta) + unitalign6[h][g][1]*cos(-zeta), unitalign6[h][g][0]*cos(-zeta) - unitalign6[h][g][1]*sin(-zeta));

                        for (n=0; n<7; n++){
                            cosk_list[n] = findClosest(cosk, cosk_single[n], 10);
                            phik_list[n] = findClosest(phik, phik_single[n], 10);
                        }
                        cosk_start = 0, cosk_end = 7;
                        phik_start = 0, phik_end = 1;
                        ang_factor = 14.3; //ie 20*5 = 100 = size(phik) * size(cosk)

                    } else {
                      continue;
                    }


                    for (n=0; n<ARRAY_SIZE; n++){  //f_polIC (compton energy)
                        for (l=0; l<ARRAY_SIZE; l++){ //f_pol (sync energy)
                            for (p=phik_start; p<phik_end; p++){ //phi
                                for (m=cosk_start; m<cosk_end; m++){ //cosk
                                    nn = m;
                                    if (h != g) {
                                       p = phik_list[nn];
                                       m = cosk_list[nn];
                                    }
                                    F_min = sqrt(f_pol_IC[n]/(2*f_pol[l] * (1-cosk[m])));
                                    v_k[0] = sqrt(1-cosk[m]*cosk[m])*cos(phik[p]), v_k[1] = sqrt(1-cosk[m]*cosk[m])*sin(phik[p]), v_k[2] = cosk[m]; //incoming photon direction vector
                                    e_k[0] = sqrt(1-cosk[m]*cosk[m])*sin(phik[p])*cos(Btheta);
                                    e_k[1] = sin(Btheta)*cosk[m] - cos(Btheta)*sqrt(1-cosk[m]*cosk[m])*cos(phik[p]); //incoming photon polarization vector (perp)
                                    e_k[2] = -sin(Btheta)*sqrt(1-cosk[m]*cosk[m])*sin(phik[p]);
                                    //this is just cross product of e_k and v_k above
                                    e_kpara[0] = pow(cosk[m],2)*sin(Btheta) - cos(Btheta)*cos(phik[p])*cosk[m]*sqrt(1-cosk[m]*cosk[m]) + (1-cosk[m]*cosk[m])*pow(sin(phik[p]),2)*sin(Btheta),
                                    e_kpara[1] = -sin(Btheta)*(1-cosk[m]*cosk[m])*sin(2*phik[p])/2 - cosk[m]*sqrt(1-cosk[m]*cosk[m])*sin(phik[p])*cos(Btheta), //incoming photon polarization vector (para)
                                    e_kpara[2] = (1-cosk[m]*cosk[m])*pow(sin(phik[p]),2)*cos(Btheta) - sqrt(1-cosk[m]*cosk[m])*cos(phik[p])*sin(Btheta)*cosk[m] + (1-cosk[m]*cosk[m])*pow(cos(phik[p]),2)*cos(Btheta);
                                    //e_k above is not normalised
                                    q_theta = pow(1-pow(cos(Btheta)*cosk[m]+sin(Btheta)*cos(phik[p])*sqrt(1-cosk[m]*cosk[m]),2),(effective_alpha[l]+1)/4);

                                    if (F_min*(Me_EV) > E_elecs[ARRAY_SIZE-1]) {
                                        P_perpIC[n] += 0.0;
                                        P_paraIC[n] += 0.0;
                                    }

                                    else if (F_min*(Me_EV) > E_elecs[0]) {
                                        for (o=0; o<ARRAY_SIZE; o++){ //find Fmin location in E_elecs
                                            if (F_min*(Me_EV) > E_elec_min[o] && F_min*(Me_EV) < E_elec_max[o]){
                                               break;
                                            }
                                        }

                                        for (i=o; i<ARRAY_SIZE; i++){ //E_e
                                            //initial constant factor to make sure this matches with isotropic electron losses, original constant from paper is 1.2859E-91

                                            Pperpperp = 3.95417e-103 /*8.26784E-118*/ * q_theta * /*dfreqs_polIC[n] */ f_pol_IC[n]/f_pol[l] * F_min * /*dfreqs_pol[l]* */  (dfactor_perp[h][g][l] /(f_pol[l]*H)) * dcosk * dphik * ang_factor //this is Power(freq), multiply by freq later to get nuF(nu)
                                            * ( Z_e(_e,e_k,v_k,1) * ( Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i])
                                            + Sigma_1(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) ) + Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) )
                                            *KN; 
                                            
                                            IC_Stokes[h][n][0] += (Pperpperp/(N_BLOCKS)), IC_Stokes[h][n][1] += (Pperpperp/(N_BLOCKS))*cos(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta))), IC_Stokes[h][n][2] += (Pperpperp/(N_BLOCKS))*sin(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta)));

                                            Pperppara = 3.95417e-103 /*8.26784E-118*/ * q_theta * /*dfreqs_polIC[n] */ f_pol_IC[n]/f_pol[l] * F_min * /*dfreqs_pol[l]* */ (dfactor_para[h][g][l] /(f_pol[l]*H)) * dcosk * dphik * ang_factor
                                            * ( Z_e(_e,e_kpara,v_k,1) * ( Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i])
                                            + Sigma_1(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) ) + Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) )
                                            *KN; 
                                            
                                            IC_Stokes[h][n][0] += (Pperppara/(N_BLOCKS)), IC_Stokes[h][n][1] += (Pperppara/(N_BLOCKS))*cos(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta))), IC_Stokes[h][n][2] += (Pperppara/(N_BLOCKS))*sin(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta)));
                                            

                                            Pparaperp = 3.95417e-103 /*8.26784E-118*/ * q_theta * /*dfreqs_polIC[n] */ f_pol_IC[n]/f_pol[l] * F_min * /*dfreqs_pol[l]* */ (dfactor_perp[h][g][l] /(f_pol[l]*H)) * dcosk * dphik * ang_factor
                                            * ( Z_e(_e,e_k,v_k,0) * ( Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i])
                                            + Sigma_1(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) ) + Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) )
                                            *KN; 
                                            
                                            IC_Stokes[h][n][0] += (Pparaperp/(N_BLOCKS)), IC_Stokes[h][n][1] += (Pparaperp/(N_BLOCKS))*cos(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta))), IC_Stokes[h][n][2] += (Pparaperp/(N_BLOCKS))*sin(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta)));

                                            

                                            Pparapara = 3.95417e-103 /*8.26784E-118*/ * q_theta * /*dfreqs_polIC[n] */ f_pol_IC[n]/f_pol[l] * F_min * /*dfreqs_pol[l]* */ (dfactor_para[h][g][l] /(f_pol[l]*H)) * dcosk * dphik * ang_factor
                                            * ( Z_e(_e,e_kpara,v_k,0) * ( Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i])
                                            + Sigma_1(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) ) + Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) )
                                            *KN; 
                                            
                                            IC_Stokes[h][n][0] += (Pparapara/(N_BLOCKS)), IC_Stokes[h][n][1] += (Pparapara/(N_BLOCKS))*cos(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta))), IC_Stokes[h][n][2] += (Pparapara/(N_BLOCKS))*sin(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta)));

                                            P_perpIC[n] += Pperpperp + Pperppara;
                                            P_paraIC[n] += Pparaperp + Pparapara;
                                            if (n == 13) {
                                                P_X[l] += Pperpperp + Pperppara + Pparaperp + Pparapara;
                                            }
                                        }
                                    }

                                    else {
                                        for (i=0; i<ARRAY_SIZE; i++){ //E_e
                                            Pperpperp = 3.95417e-103 /*8.26784E-118*/ * q_theta * /*dfreqs_polIC[n] */ f_pol_IC[n]/f_pol[l] * F_min * /*dfreqs_pol[l]* */  (dfactor_perp[h][g][l] /(f_pol[l]*H)) * dcosk * dphik * ang_factor
                                            * ( Z_e(_e,e_k,v_k,1) * ( Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i])
                                            + Sigma_1(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) ) + Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) )
                                            *KN; 
                                            
                                            IC_Stokes[h][n][0] += (Pperpperp/(N_BLOCKS)), IC_Stokes[h][n][1] += (Pperpperp/(N_BLOCKS))*cos(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta))), IC_Stokes[h][n][2] += (Pperpperp/(N_BLOCKS))*sin(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta)));

                                            Pperppara = 3.95417e-103 /*8.26784E-118*/ * q_theta * /*dfreqs_polIC[n] */ f_pol_IC[n]/f_pol[l] * F_min * /*dfreqs_pol[l]* */ (dfactor_para[h][g][l] /(f_pol[l]*H)) * dcosk * dphik * ang_factor
                                            * ( Z_e(_e,e_kpara,v_k,1) * ( Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i])
                                            + Sigma_1(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) ) + Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) )
                                            *KN; 
                                            
                                            IC_Stokes[h][n][0] += (Pperppara/(N_BLOCKS)), IC_Stokes[h][n][1] += (Pperppara/(N_BLOCKS))*cos(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta))), IC_Stokes[h][n][2] += (Pperppara/(N_BLOCKS))*sin(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta)));
                                            

                                            Pparaperp = 3.95417e-103 /*8.26784E-118*/ * q_theta * /*dfreqs_polIC[n] */ f_pol_IC[n]/f_pol[l] * F_min * /*dfreqs_pol[l]* */ (dfactor_perp[h][g][l] /(f_pol[l]*H)) * dcosk * dphik * ang_factor
                                            * ( Z_e(_e,e_k,v_k,0) * ( Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i])
                                            + Sigma_1(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) ) + Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) )
                                            *KN; 
                                            
                                            IC_Stokes[h][n][0] += (Pparaperp/(N_BLOCKS)), IC_Stokes[h][n][1] += (Pparaperp/(N_BLOCKS))*cos(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta))), IC_Stokes[h][n][2] += (Pparaperp/(N_BLOCKS))*sin(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta)));

                                            

                                            Pparapara = 3.95417e-103 /*8.26784E-118*/ * q_theta * /*dfreqs_polIC[n] */ f_pol_IC[n]/f_pol[l] * F_min * /*dfreqs_pol[l]* */ (dfactor_para[h][g][l] /(f_pol[l]*H)) * dcosk * dphik * ang_factor
                                            * ( Z_e(_e,e_kpara,v_k,0) * ( Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i])
                                            + Sigma_1(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) ) + Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) )
                                            *KN; 
                                            
                                            IC_Stokes[h][n][0] += (Pparapara/(N_BLOCKS)), IC_Stokes[h][n][1] += (Pparapara/(N_BLOCKS))*cos(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta))), IC_Stokes[h][n][2] += (Pparapara/(N_BLOCKS))*sin(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta)));

                                            P_perpIC[n] += Pperpperp + Pperppara;
                                            P_paraIC[n] += Pparaperp + Pparapara;
                                            if (n == 13) {
                                                P_X[l] += Pperpperp + Pperppara + Pparaperp + Pparapara;
                                            }

                                        }
                                    }

                                    if (h != g) {
                                       m = nn;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for (i=0; i<ARRAY_SIZE; i++){
            for (n=0; n<N_BLOCKS; n++)
                u_rad += (urad_array_perp[n][i]+urad_array_para[n][i]);
        }
        u_rad /= N_BLOCKS; //use average energy density for electron losses, otherwise different blocks cool at different rates - SSC cooling isnt that important anyway
        
        for (i=0; i<ARRAY_SIZE; i++) //1.79 //0.85
        {
            Ps_per_mIC_elec[i] = 0.0;
            Ps_per_mIC_elec[i] = (4/3)*6.6524E-29*gamma_e[i]*gamma_e[i]*beta_e[i]*beta_e[i]*dN_dE[i]*dEe[i]*Qe*u_rad; //test for Ps_per_mIC_elec == this
            IC_losses[i] = Ps_per_mIC_elec[i];
        }
        
        
        //**********************************************************************************************//  
        //some code to estimate what dx should be
        dx_set = findminelementNO0(Ne_e, ARRAY_SIZE);//gets the lowest non_zero element. Higher energy electrons radiate more rapidly


        dx_R = 0.05*(R0+x*tan(deg2rad(theta_open_p)))/tan(deg2rad(theta_open_p)); //ensures Rnew <= 1.05 Rold
        dx_P = (Ne_e[dx_set]*(E_elecs[dx_set]-E_elecs[dx_set-1])*Qe)/(Ps_per_m_elec[dx_set] + Ps_per_mIC_elec[dx_set]); //based on radiative losses
  
        if (dx_R < dx_P) {
            dx = dx_R;
        } else {
            dx = dx_P;
        }

        for (h=0; h<N_BLOCKS; h++){ //summing the Stokes parameters of all the blocks and boosting
            zone_doppler = doppler(beta_bulk, theta_tot[h]);
            //difference between bin boundaries in logspace:
            logdif = log10(f_pol_IC[1])-log10(f_pol_IC[0]);
            binshiftIC = (int)round(log10(zone_doppler/doppler_factor)/logdif); //the number of bins one moves across (relative to shift of bins which already takes place)
            logdif = log10(f_pol[1])-log10(f_pol[0]);
            binshiftS = (int)round(log10(zone_doppler/doppler_factor)/logdif);
            printf("binshifts %d\t%d\n", binshiftIC,binshiftS);
            printf("dopplers %.5e\t%.5e\t%.5e\n", zone_doppler,doppler_factor,logdif);

            for (n=0; n<ARRAY_SIZE; n++){
                // Might want zone_doppler ^ 3 for a continuous jet
                if (n-binshiftIC >= 0 && n-binshiftIC < ARRAY_SIZE){
                    IC_StokesTotal[n][0] += IC_Stokes[h][n-binshiftIC][0] * pow(zone_doppler,4) * dx;
                    IC_StokesTotal[n][1] += IC_Stokes[h][n-binshiftIC][1] * pow(zone_doppler,4) * dx;
                    IC_StokesTotal[n][2] += IC_Stokes[h][n-binshiftIC][2] * pow(zone_doppler,4) * dx;

                    IC_StokesZTotal[h][n][0] += IC_Stokes[h][n-binshiftIC][0] * pow(zone_doppler,4) * dx;
                    IC_StokesZTotal[h][n][1] += IC_Stokes[h][n-binshiftIC][1] * pow(zone_doppler,4) * dx;
                    IC_StokesZTotal[h][n][2] += IC_Stokes[h][n-binshiftIC][2] * pow(zone_doppler,4) * dx;
                    
                }
                if (n-binshiftS >= 0 && n-binshiftS < ARRAY_SIZE){
                    S_StokesTotal[n][0] += S_Stokes[h][n-binshiftS][0] * pow(zone_doppler,4) * dx; 
                    S_StokesTotal[n][1] += S_Stokes[h][n-binshiftS][1] * pow(zone_doppler,4) * dx; 
                    S_StokesTotal[n][2] += S_Stokes[h][n-binshiftS][2] * pow(zone_doppler,4) * dx; 

                    S_StokesZTotal[h][n][0] += S_Stokes[h][n-binshiftS][0] * pow(zone_doppler,4) * dx; 
                    S_StokesZTotal[h][n][1] += S_Stokes[h][n-binshiftS][1] * pow(zone_doppler,4) * dx; 
                    S_StokesZTotal[h][n][2] += S_Stokes[h][n-binshiftS][2] * pow(zone_doppler,4) * dx; 
                }
            }
        }

        //change values for next loop
        x += dx;
        R_prev = R;
        R = R_new(R0, deg2rad(theta_open_p), x);

        //correct synchrotron for full length + IC
        for (i=0; i<ARRAY_SIZE; i++)
        {
            Ps_per_m[i] *= dx;
            //Ps_per_m_IC[i] *= dx;
            Sync_losses[i] *= dx;
            IC_losses[i] *= dx;              //ALP
            P_perp[i] *= dx;  //polarisation powers
            P_para[i] *= dx;
            P_perpIC[i] *= dx;
            P_paraIC[i] *= dx;
            P_X[i] *= dx;
	    //fprintf(Xfile,"\t%.5e", P_X[i] * f_pol[i]);
            //fprintf(Proj_Bfile,"\t%.5e",Proj_theta_B);

            if (i==ARRAY_SIZE-1) //puts the outputs for the next jet section on a new line in files
            {  
               //fprintf(Xfile,"\n");
               //fprintf(Proj_Bfile, "\n");
            }
        }

        //compute the synchrotron losses + IC LOSSES
        for (i=0; i<ARRAY_SIZE; i++) {
            if (i==0) {
                //Ne_losses[i]=0.0; //cannot drop to lower bins
                Ne_losses[i]=(Sync_losses[i] + IC_losses[i]) / ((E_elecs[i]-Me_EV)*Qe); //can't radiate anymore
            }
            else {
	      //Ne_losses[i] = (Sync_losses[i])*(C/C)/((E_elecs[i]-E_elecs[i-1])*Qe);
	      Ne_losses[i] = (Sync_losses[i] + IC_losses[i]) / ((E_elecs[i]-E_elecs[i-1])*Qe);
            }
            if (i==ARRAY_SIZE-1) {
                Ne_gains[i] = 0.0;
            }
            else {
	      //Ne_gains[i] = (Sync_losses[i+1])*(C/C)/((E_elecs[i+1]-E_elecs[i])*Qe);
	      Ne_gains[i] = (Sync_losses[i+1] + IC_losses[i+1])*(C/C)/((E_elecs[i+1]-E_elecs[i])*Qe);
            }
        }

        for (i=0; i<ARRAY_SIZE; i++) {
            Ne_e[i] = Ne_e[i] - Ne_losses[i] + Ne_gains[i];
            if (x>4.9E12 && intermed_param==0) {
                Ne_intermed[i] = Ne_e[i];
                if (i==ARRAY_SIZE-1) {
                    intermed_param = 1; //ensures condition only met once
                }
            }

        }

        cleanpop(Ne_e, ARRAY_SIZE, 10.0); //any elements with Ne_e<10 set to zero to avoid nans
        cleanpop(Ne_intermed, ARRAY_SIZE, 10.0); //any elements with Ne_e<10 set to zero to avoid nans
                
	B_prev = B;//need to update dfreqs
        B = get_newB(B0, R0, R);
        eps = epsilon(B);

        // reset values to 0 for the next section x step
	u_rad = 0.0;
        for (i=0; i<ARRAY_SIZE; i++)
        {
            f_c[i] = f_crit(gamma_e[i], B);
            f_em_min[i] *= B/B_prev;//steps[nSteps];
            f_em_max[i] *= B/B_prev;//steps[nSteps];
            dfreqs[i] = f_em_max[i]-f_em_min[i];
            A_elecs[i] = A_orig*(Ne_e[i]/Ne_orig[i]);
            dN_dE[i] = electron_PL(A_elecs[i], alpha, E_elecs[i]*Qe, E_max*Qe);
            P_para[i] = 0.0;//Reset
            P_perp[i] = 0.0;
            P_paraIC[i] = 0.0;
            P_perpIC[i] = 0.0;
            P_X[i] = 0.0;
            Ps_per_mIC_elec[i] = 0.0;
            Ps_per_m_elec[i] = 0.0;

        }

        memset(S_Stokes, 0, sizeof(S_Stokes[0][0][0])* N_BLOCKS * ARRAY_SIZE * 3); //reset these to 0 for next section
        memset(IC_Stokes, 0, sizeof(IC_Stokes[0][0][0])* N_BLOCKS * ARRAY_SIZE * 3);
        memset(urad_array_perp, 0, sizeof(urad_array_perp[0][0])* N_BLOCKS * ARRAY_SIZE); //reset these to 0 for next section
        memset(urad_array_para, 0, sizeof(urad_array_para[0][0])* N_BLOCKS * ARRAY_SIZE);
        memset(dfactor_perp, 0, sizeof(dfactor_perp[0][0][0])* N_BLOCKS * N_BLOCKS * ARRAY_SIZE); //reset these to 0 for next section
        memset(dfactor_para, 0, sizeof(dfactor_para[0][0][0])* N_BLOCKS * N_BLOCKS * ARRAY_SIZE);

        nSteps+=1; //allows the number of jet sections to be determined
        printf("nStep, dx, x, B, R %.5d\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n", nSteps, dx, x, B, R, S_StokesTotal[14][0] * f_pol[14], IC_StokesTotal[13][0] * f_pol_IC[13]);

        //save data needed at every jet section to solve line of sight opacity
        fprintf(basicdata, "\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e \n", dx, x, B, R, S_StokesTotal[14][0] * f_pol[14], IC_StokesTotal[13][0] * f_pol_IC[13]);
        
        /************************************/
        //Timing: to understand accurate cluster submission
        time_end = clock();
        time_spent = (double)(time_end - time_begin) / CLOCKS_PER_SEC;
        printf("Time spent this section: %.3f mins \n", time_spent/60);
        /************************************/

        if (nSteps == breakstep) {
            break;
        }
    } // End of Jet calculation Loop ----------------------------------------------------------------------------------------//


    free(R_array);
    free(dx_array); //free all dynamically allocated memory
    for (i=0; i<buffer_size; i++) {
	free(Pperp_array[i]);
        free(Ppara_array[i]);
    }
    free(Pperp_array);
    free(Ppara_array);

    double newbins[ARRAY_SIZE]; // rebinning for IC to fit into synchrotron bins, choose closest synchrotron flux to be added to IC flux
    for (n=0; n<ARRAY_SIZE; n++){
        ICS_StokesTotal[n][0] = S_StokesTotal[n][0] + IC_StokesTotal[findClosest(f_pol_IC, f_pol[n], ARRAY_SIZE)][0];
        ICS_StokesTotal[n][1] = S_StokesTotal[n][1] + IC_StokesTotal[findClosest(f_pol_IC, f_pol[n], ARRAY_SIZE)][1]; // ICS uses IC binning
        ICS_StokesTotal[n][2] = S_StokesTotal[n][2] + IC_StokesTotal[findClosest(f_pol_IC, f_pol[n], ARRAY_SIZE)][2];
    }

    for (n=0; n<ARRAY_SIZE; n++){ //multiplying F(nu) by nu to get nuF(nu), to plot Power instead: switch f_pol by dfreqs_pol, however this will require a different treatment of ICS_total
        ICS_StokesTotal[n][0] *= f_pol[n];
        ICS_StokesTotal[n][1] *= f_pol[n];
        ICS_StokesTotal[n][2] *= f_pol[n];

        IC_StokesTotal[n][0] *= f_pol_IC[n];
        IC_StokesTotal[n][1] *= f_pol_IC[n];
        IC_StokesTotal[n][2] *= f_pol_IC[n];

        S_StokesTotal[n][0] *= f_pol[n];
        S_StokesTotal[n][1] *= f_pol[n];
        S_StokesTotal[n][2] *= f_pol[n];

        for (h=0; h<N_BLOCKS; h++){
             IC_StokesZTotal[h][n][0] *= f_pol_IC[n];
             IC_StokesZTotal[h][n][1] *= f_pol_IC[n];
             IC_StokesZTotal[h][n][2] *= f_pol_IC[n];

             S_StokesZTotal[h][n][0] *= f_pol[n];
             S_StokesZTotal[h][n][1] *= f_pol[n];
             S_StokesZTotal[h][n][2] *= f_pol[n];
        }
    }

    for (i = 0; i < N_BLOCKS; i++)
    {
      for (m = 0; m < ARRAY_SIZE; m++)
      {
        for (n = 0; n < 3; n++)
        {
          fprintf(IC_Z, "%f ", IC_StokesZTotal[i][m][n]);
          fprintf(S_Z, "%f ", S_StokesZTotal[i][m][n]);
        }
        fprintf(IC_Z, "\n");
        fprintf(S_Z, "\n");
      }
      fprintf(IC_Z, "\n");
      fprintf(S_Z, "\n");
    }


    for (n=0; n<ARRAY_SIZE; n++){
            IC_Pi[n] = sqrt(IC_StokesTotal[n][1]*IC_StokesTotal[n][1] + IC_StokesTotal[n][2]*IC_StokesTotal[n][2]) / IC_StokesTotal[n][0];
            S_Pi[n] = sqrt(S_StokesTotal[n][1]*S_StokesTotal[n][1] + S_StokesTotal[n][2]*S_StokesTotal[n][2]) / S_StokesTotal[n][0];
            ICS_Pi[n] = sqrt(ICS_StokesTotal[n][1]*ICS_StokesTotal[n][1] + ICS_StokesTotal[n][2]*ICS_StokesTotal[n][2]) / ICS_StokesTotal[n][0];

            IC_PA[n] = 0.5*atan2(IC_StokesTotal[n][2],IC_StokesTotal[n][1]);
            S_PA[n] = 0.5*atan2(S_StokesTotal[n][2],S_StokesTotal[n][1]);
            ICS_PA[n] = 0.5*atan2(ICS_StokesTotal[n][2],ICS_StokesTotal[n][1]);

            IC_P[n] = IC_StokesTotal[n][0];
            S_P[n] = S_StokesTotal[n][0];
            ICS_P[n] = ICS_StokesTotal[n][0];

            fprintf(PI, "\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n",S_Pi[n],S_PA[n],
                   S_P[n],IC_Pi[n],IC_PA[n],IC_P[n],ICS_Pi[n],ICS_PA[n],ICS_P[n]);

    }

    printf("The jet was divided into %d %s \n", nSteps, "section(s).");
    return 0;
}


