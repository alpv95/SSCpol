#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "jet_fns.h"
#include "mtwister.h"
#include <time.h>



/* Program generating quiescent synchrotron and IC (thomson limit) power spectrum and polarization for thin jet
* Includes Polarization of all components
*
*
*
*
*
*
*/

//define some physical constants
const double Me_EV=0.511E6; //electron rest energy
const double C = 3.0E8; //speed of light
const double Qe = 1.6E-19; //electron charge
const double Me = 9.11E-31; //electron mass kg
const double H = 6.63E-34; //Plancks constant



//define some looping parameters
int array_size=50; // sets the number of synchrotron & IC pts
double array_size_d=50.0; // use to set log ratios, should be the same as above line but with .0
int i, l, m, n, o, p, g, h; //some looping parameters
double w, x, y, z, a, b, c, d, q;
int dx_set; // use to define smallest non zero population

int main(int argc,char* argv[]) //argc is integer number of arguments passed, argv[0] is program name, argv[1..n] are arguments passed in string format
{

    //enter the jet parameters
    double W_j = 3E37;  //3e37//W_jarray[counter]; //W jet power in lab frame 7.63E36. Should be OBSERVED POWER
     double L_jet = 5E20;// 6E11;//1E19;//6E12;// 1E19;// 6E12;//6E12; //1.0E19; //length in m in the fluid frame
     double B = 4E-5, B0 = 4E-5; //4e-3 4e-5 //8e-5//0.000268; //B-field at jet base
     double R0 = 0.0, R = 0.0; //4.532E13;//7.32E13; // Radius of the jet at the base 3.32 works fairly well 
    double R_prev = 0.0, B_prev = 0.0; //changing parameters of the jet-initialise. R prev corrects for increasing jet volume 
    double E_min = 5.11E6; // Minimum electron energy 
    double E_max = 8.1E9;//2.5E10 8.1E9;//50e9//8.1E9//5.0E9;//5.60E9; // Energy of the ECO in eV 
    double alpha = 1.95;//1.95;//1.9//1.95//2.000001; // PL index of electrons
     double theta_open_p = 45.0;// 60 50//*(M_PI/180.0); // opening angle of the jet in the fluid frame 
    double theta_obs; //4//3.0;//*(M_PI/180.0); // observers angle to jet axis in rad 
    sscanf(argv[4], "%lf", &theta_obs);
    double gamma_bulk = 7.5;//pow(10,(log10(W_j)*0.246-8.18765 + 0.09)); //final additive constant to make sure highest is 40 and lowest is 5//12.0; // bulk Lorentz factor of jet material
    int n_blocks;//127; //for the TEMZ model, can have 1,7,19,37,61,91,127 blocks, (rings 0,1,2,3,4,5,6)
    sscanf(argv[5], "%d", &n_blocks);
    int n_rings; //6; //up to 6 rings possible atm, must choose number of rings corresponding to number of zones
    sscanf(argv[6], "%d", &n_rings);
    printf("%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%d\n", W_j, gamma_bulk, theta_obs, B, B0, E_max, alpha, n_blocks);



    //define radius at jet base
    R0 = R_0(W_j*(3.0/4.0), 1.0, gamma_bulk, B0);//7.32E13;
    R = R_0(W_j*(3.0/4.0), 1.0, gamma_bulk, B0);//7.32E13;
    printf("Radius at jet base: %.5e \n", R);

    //read in any required data (only Bessel Functions for now)
    FILE *FGfile;
    FGfile = fopen("FG.txt","r"); //Bessel Functions for Synchrotron

    double FG[75][3];
    for (m=0; m<75; m++)
    {
           fscanf(FGfile, "%lf\t%lf\t%lf\n", &q, &c, &d);
           FG[m][0] = q;
           FG[m][1] = c;
           FG[m][2] = d;
    }

    //define some files to store output data
    FILE *Nfile, *freqrange, *freqfile, *opfile, *powfile, *basicdata, *elecpop, *elecpopfull, *photonpop, *jfile, *keyparams, *ICfile,
        *Pperpfile, *Pparafile, *PperpfileIC, *PparafileIC, *Proj_Bfile, *block_thetafile, *TESTFIL2;

    Nfile = fopen("Ndata.txt", "w"); //store electron populations
    freqrange = fopen("freqrange.txt", "w");//store frequency bin boundaries
    freqfile = fopen("critfreqs.txt", "w");//store critical frequencies
    opfile = fopen("kvalues.txt", "w");//store opacity values for each jet section
    powfile = fopen("secpow.txt", "w");//store emitted synchrotron power (freq fn) for each jet section
    basicdata = fopen("basicdata.txt", "w"); //store B, x, R etc
    elecpop = fopen("elecpop.txt", "w"); //store electron populations at various points down the jet
    elecpopfull = fopen("elecpopfull.txt", "w"); //store electron populations at all points down the jet
    photonpop = fopen("photonpop.txt", "w"); //store sync photon populations at various points down the jet
    jfile = fopen("jfile_re.txt", "w");
    keyparams = fopen("keyparams.txt", "w"); //store Lj, gamma_bulk, theta_obs and theta_open
    ICfile = fopen("ICfile.txt", "w"); //ALP
    Pperpfile = fopen("Pperpfile.txt", "w"); //polarisations
    Pparafile = fopen("Pparafile.txt", "w");
    PperpfileIC = fopen("PperpfileIC.txt", "w"); //polarisations
    PparafileIC = fopen("PparafileIC.txt", "w");
    Proj_Bfile = fopen("Proj_Bfile.txt", "w"); //projected B field onto plane of the sky for each section
    block_thetafile = fopen("block_thetafile.txt","w"); //saves the angle to the line of sight of each block
    TESTFIL2 = fopen("TESTFIL2.txt","w");
    fprintf(keyparams, "\t%.5e\t%.5e\t%.5e\t%.5e \n", L_jet, gamma_bulk, theta_obs, theta_open_p);

    //define some useful parameters to gauge progress
    int nSteps = 0; //how many sections the jet has broken down into
    double x = 0; //progress along the jet
    double dx; //incremenet of step length
    double eps; //stores the eqn 2.18b
    double beta_bulk = lortovel(gamma_bulk);
    double doppler_factor = doppler(beta_bulk, deg2rad(theta_obs));
    int intermed_param = 0; //allows intermediate population to be determined to compare to paper

    //define some parameters for determining the dx of each jet section
    double dx_R; //dx where R_new = 1.05*R_old
    double dx_P; //dx to ensure smallest bin population not depleted

    //some parameters used to define the population of electrons
    double A_elecs[array_size]; //stores PL coefficient
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
    double Ee_max= 4.0E12;//Me_EV*3.0E4;//4.0E12; //sets the maximum electron energy in eV
    double ratio = pow((Ee_max/E_min), (1/array_size_d)); //logarithmic ratio of electron bins: USED FOR BOUNDS
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
    double dfreqs[array_size];// store the bin widths in critical frequency
    double f_em_min[array_size], f_em_max[array_size], f_em_min_IC[array_size], f_em_max_IC[array_size]; // store bin boundaries for critical frequencies
    double j[array_size]; // emissivity per Hz per unit volume
    double k[array_size]; // opacity, units of m^-1
    double Ps_per_m[array_size]; // Synchrotron power per unit m
    double Ps_per_m_elec[array_size]; //power lost by each electron bin per m
    double Ps_per_m_test[array_size]; //used to normalise IC power emission
    memset(Ps_per_m_test, 0.0, array_size*sizeof(Ps_per_m_test[0]));
    double Sync_losses[array_size];// store total electron energy losses
    double P_single[array_size];// power emitted by a single electron of energy E
    double Ps_per_m_IC[array_size]; //IC power per m
    double Ps_per_mIC_elec[array_size]; //energy lost by each electron bin per m of IC emission
    memset(Ps_per_mIC_elec, 0.0, array_size*sizeof(Ps_per_mIC_elec[0]));
    double urad_array_perp[array_size]; //photon energy density
    double urad_array_para[array_size];
    double IC_losses[array_size]; //store electron IC energy losses

    //****************Initialise Polarisation variables *****************************************************************************//
    //first set frequency bins for polarised powers to drop into:
    double f_pol[array_size];
    double f_pol_IC[array_size]; //fixed bins
    double f_pol_ICbounds[array_size+1]; //for calculating IC bins properly
    double P_perp[array_size];
    memset(P_perp, 0.0, array_size*sizeof(P_perp[0]));
    double P_para[array_size];
    memset(P_para, 0.0, array_size*sizeof(P_para[0])); //initialise to make sure junk values arent inside upon first +=
    double P_perpIC[array_size];
    memset(P_perpIC, 0.0, array_size*sizeof(P_perpIC[0]));
    double P_paraIC[array_size];
    memset(P_paraIC, 0.0, array_size*sizeof(P_paraIC[0]));
    double dfreqs_pol[array_size]; //again fixed size
    double dfreqs_polIC[array_size];
    double F_min;
    //cos(th_k) for compton integral
    double cosk[10];
    double phik[10];
    double dcosk = 1/5.5; //1/25.5;
    double dphik = (2*M_PI)/10;//(2*M_PI)/array_size;
    for (i=0; i<10; i++){
        cosk[i] = (i-5)/5.5; //(i-25)/25.5;
        phik[i] = -M_PI + i*(2*M_PI)/10; //i*(2*M_PI)/array_size;
    }


    //Fill electron arrays with initial parameters
    for (i=0; i<array_size; i++)
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

    }//ends electron population definition.

    //sets equivalent IC bins and saves freqrange
    //log10spaceWithMidpoints(f_pol_ICbounds, f_pol_IC, gamma_e[0] * gamma_e[0] * f_pol[0],
                            //gamma_e[array_size-1] * gamma_e[array_size-1] * f_pol[array_size-1], array_size+1); //sets IC frequency bins and bounds given start and end points
    double bindiff;
    for (i=0; i<array_size; i++) {
        bindiff = log10(gamma_e[1] * gamma_e[1] * (f_em_min[1]))-log10(gamma_e[0] * gamma_e[0] * (f_em_max[0]));
        f_em_min_IC[i] = pow(10,log10(gamma_e[i] * gamma_e[i] * f_em_min[i])-bindiff/2);//pow(10,(log10(f_pol_ICbounds[i])+0.15));
        f_em_max_IC[i] = pow(10,log10(gamma_e[i] * gamma_e[i] * f_em_max[i])+bindiff/2);//pow(10,(log10(f_pol_ICbounds[i+1])-0.15));
        dfreqs_polIC[i] = (f_em_max_IC[i]-f_em_min_IC[i]);
        f_pol_IC[i] = pow(10,log10(gamma_e[i] * gamma_e[i] * f_pol[i]));

        fprintf(freqrange,  "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e \n", f_em_min[i], f_em_max[i], f_c[i], dfreqs[i],
                f_em_min_IC[i], f_em_max_IC[i], f_pol_IC[i], dfreqs_polIC[i]);  //saves the frequency bins!

    }


    A_orig = A_elecs[0]; //all elements the same to begin with


    //print the energies to test equipartition
    UB_jetbase = (4.0/3.0)*M_PI*R0*R0*C*B0*B0/(2*4*M_PI*1E-7);
    //printf("Obs power, particles, B-fields: %.5e\t%.5e\t%.5e\t%.5e \n", W_j/(gamma_bulk*gamma_bulk), Ue_jetbase, UB_jetbase, Ue_jetbase+UB_jetbase);
    
    //redefine electron population to ensure equipartition
    for (i=0; i<array_size; i++)
    {

        A_elecs[i] = A_PL(alpha, W_j, gamma_bulk, E_min*Qe, E_max*Qe, 1.0)*UB_jetbase/Ue_jetbase; //1 assumes equipartition. DIVIDED W_j by gamma^2 for jet frame
        dN_dE[i] = electron_PL(A_elecs[i], alpha, E_elecs[i]*Qe, E_max*Qe);
        Ne_e[i] = electron_PL(A_elecs[i], alpha, E_elecs[i]*Qe, E_max*Qe)*dEe[i]*Qe;
        Ne_orig[i] = electron_PL(A_elecs[i], alpha, E_elecs[i]*Qe, E_max*Qe)*dEe[i]*Qe; //allows comparison of final to initial electron population

    }//ends electron population re-definition.
    
    //recalaculate equipartition fraction
    Ue_jetbase = 0.0; 
    for (i=0; i<array_size; i++)
    {
        Ue_jetbase += Ne_e[i]*E_elecs[i]*Qe;
    }

    //fudge factors: set bins with fewer than ten electrons to zero: prevents resolution from being limited
    cleanpop(Ne_e, array_size, 10.0);
    cleanpop(dN_dE, array_size, 10.0);
    

    for (i=0; i<array_size; i++)
    {
        A_elecs[i]*=(Ne_e[i]/Ne_orig[i]); //sets A to zero for tiny populations
    }

    
    //11/10/16 define some parameters to determine thick/thin jet sections
    double R_eff[array_size]; //effective radius down to which can be seen -> smaller than jet radius if optically thick


    //HELICAL B-field: This variable will save the 3D B-fields on the sky for each jet zone, these B-fields evolve along the jet
    //becoming more transverse
    double Proj_theta_B[n_blocks];
    //defining how far along the helix we are
    int helix_counter;
    sscanf(argv[2], "%d", &helix_counter); //taking input from command line how many 'days' we have been observing
    //parameters for helix
    double t_helix = helix_counter*M_PI/32; //helix parameter x=rcos(t+phi) y=rsin(t+phi) z =ct (r will be radius of jet)
    //Pi constant decided how quickly we sample along the helix -> timescale
    double phi_helix = 0;//M_PI/3; //phase shift, just an initial condition
    double c_helix = R0; //speed (ie number of coils per unit distance) -- high c => spaced out coils, low c => densely packed coils

    //choosing random vectors (B_0,B_1,B_2) in unit sphere for random blocks initial B directions:
    MTRand seedr = seedRand((unsigned)time(NULL));
    double B_0[n_blocks];
    double B_1[n_blocks];
    double B_2[n_blocks];
    double B_length;
    for (i=0; i<(n_blocks); i++){ //these are the random B-field vectors in each block as if we are looking straight down jet (have to rotate by theta_obs)
      B_0[i] = genRand(&seedr)*2;
      B_1[i] = genRand(&seedr)*2;
      B_2[i] = genRand(&seedr)*2;
      //printf("a, b %.5e\t%.5e\n", B_0[i], B_1[i]);
      B_0[i]=(B_0[i]-1);
      B_1[i]=(B_1[i]-1);
      B_2[i]=(B_2[i]-1);
      B_length = sqrt(B_0[i]*B_0[i]+B_1[i]*B_1[i]+B_2[i]*B_2[i]);
      B_0[i] = B_0[i] / (B_length); //Normalising B-field vectors
      B_1[i] = B_1[i] / (B_length);
      B_2[i] = B_2[i] / (B_length);
      //printf("randBs %.5e\t%.5e\t%.5e\n",B_0[i],B_1[i],B_2[i]);
    }


    int EVPA_rotation;
    sscanf(argv[1], "%d", &EVPA_rotation); //taking input from command line if true begin EVPA rotation, if false do not
    int DopDep;
    sscanf(argv[3], "%d", &DopDep); //taking input from command line if true activate DopDep, if false do not

    //now the Bfield vectors are rotated about both the y and x axis each block with a different theta_obs
    //these give the angles each block is at with respect to the direction of the jet:
    //theta_x => rotation angle about x axis, theta_y => rotation angle about y axis
    double th = 2*atan(tan(deg2rad(theta_open_p))/gamma_bulk); //full jet diameter opening angle in lab frame in rad
    double theta_r[n_blocks];                // Marscher configuration, rings of turbulent cells, middle 3 rings turbulent
    double theta_phi[n_blocks];
    //TEMZ config, calculates theta_r and theta_phi for each block in the concentric circles
    for (i=0; i<n_blocks; i++){
        if (i<1){
            theta_r[i] = 0;
            theta_phi[i] = 0;
        }
        else if (i>0 && i<7){

            theta_r[i] = th/(2*n_rings);
            theta_phi[i] = -M_PI+(i-1)*2*M_PI/6;

        }
        else if (i>6 && i<19){

            theta_r[i] = 2*th/(2*n_rings);
            theta_phi[i] = -M_PI+(i-7)*2*M_PI/12;

        }
        else if (i>18 && i<37){

            theta_r[i] = 3*th/(2*n_rings);
            theta_phi[i] = -M_PI+(i-19)*2*M_PI/18;

        }
        else if (i>36 && i<61){

            theta_r[i] = 4*th/(2*n_rings);
            theta_phi[i] = -M_PI+(i-37)*2*M_PI/24;

        }
        else if (i>60 && i<91){

            theta_r[i] = 5*th/(2*n_rings);
            theta_phi[i] = -M_PI+(i-61)*2*M_PI/30;

        }
        else if (i>90 && i<127){

            theta_r[i] = 6*th/(2*n_rings);
            theta_phi[i] = -M_PI+(i-91)*2*M_PI/36;

        }
    }

    double theta_tot[n_blocks]; //total angle to line of sight of each block (for doppler boosting) in radians
    for (i=0; i<n_blocks; i++){
        if (n_blocks==1){
            theta_tot[i] = deg2rad(theta_obs);
        } else {
            theta_tot[i] = acos(cos(theta_r[i])*cos(deg2rad(theta_obs)) + sin(theta_r[i])*sin(deg2rad(theta_obs))*cos(M_PI - fabs(theta_phi[i])));
        }
    }

    for (i=0; i<n_blocks; i++){
        fprintf(block_thetafile, "\t%.5e", theta_tot[i]); //saving the angle to the direction of view of each block to be used in python file for doppler boost
    }
    //alignment vector matrix for each block to every other block, for use in Synchrotron mixing between zones for IC emission:
    //assuming jet section is flat 2D (ok assumption given large section size) but have z component to rotate with theta_obs
    //remember to add +1 to vector lengths so that 1/r is never 1/0
    double align[n_blocks][n_blocks][3]; //vector from m zone to i zone
    for (i=0; i<n_blocks; i++){
        for (m=0; m<n_blocks; m++){
            align[i][m][0] = cos(deg2rad(theta_obs))*((theta_r[i]/(th/(2*n_rings)))*cos(theta_phi[i]) - (theta_r[m]/(th/(2*n_rings)))*cos(theta_phi[m]));
            align[i][m][1] = (theta_r[i]/(th/(2*n_rings)))*sin(theta_phi[i]) - (theta_r[m]/(th/(2*n_rings)))*sin(theta_phi[m]);
            align[i][m][2] = -sin(deg2rad(theta_obs))*((theta_r[i]/(th/(2*n_rings)))*cos(theta_phi[i]) - (theta_r[m]/(th/(2*n_rings)))*cos(theta_phi[m]));
        }
    }


    double Bx[n_blocks];
    double By[n_blocks];
    double Bz[n_blocks];
    double Bx_helical[n_blocks]; //simplifies angular shifts for the helical components
    double By_helical[n_blocks];
    double Bz_helical[n_blocks];
    double BX[n_blocks]; //for joint helical and random array
    double BY[n_blocks];
    double BZ[n_blocks];

    //Sectional, or possibly full, Stokes vector bins for Synchrotron and IC respectively
    double S_Stokes[n_blocks][array_size][3];
    memset(S_Stokes, 0, sizeof(S_Stokes[0][0][0])* n_blocks * array_size * 3);
    double S_StokesTotal[array_size][3];
    memset(S_StokesTotal, 0, sizeof(S_StokesTotal[0][0])* array_size * 3);
    double S_Pi[array_size];
    double S_PA[array_size];
    double S_P[array_size];

    double IC_Stokes[n_blocks][array_size][3];
    memset(IC_Stokes, 0, sizeof(IC_Stokes[0][0][0])* n_blocks * array_size * 3);
    double IC_StokesTotal[array_size][3];
    memset(IC_StokesTotal, 0, sizeof(IC_StokesTotal[0][0])* array_size * 3);
    double IC_Pi[array_size];
    double IC_PA[array_size];
    double IC_P[array_size];
    double zone_doppler;
    double q_theta;
    double zeta;
    double psi;
    int phik_start, phik_end, cosk_start, cosk_end;
    double cosk_single, phik_single;
    double u_rad = 0.0;
    double Blength; //for normalization
    double logdif;
    int binshiftIC;
    int binshiftS;
    double dfactor;

    double X = 0.0;
    double KN = 1.0;
    double v_k[3] = {0};
    double e_k[3] = {0};
    double _e[2] = {0};
    double e_kpara[3] = {0};
    double zone_d; //distance from other zones
    //now including full RPAR treatment for IC
    double Btheta = 90*M_PI/180;
    double Pperpperp = 0.0;
    double Pperppara = 0.0;
    double Pparaperp = 0.0;
    double Pparapara = 0.0;
    double effective_alpha[array_size];



    //B_effectives due to RPAR for each zone in every other zone
    double B_effectives[n_blocks][n_blocks][3];
    memset(B_effectives, 0, sizeof(B_effectives[0][0][0])* n_blocks * n_blocks * 3);

    printf("Beginning jet analysis... \n");
    while (x<L_jet)//for (x=0; x<L_jet; x+=dx)
    {
        //update jet parameters
        eps = epsilon(B);//same for all populations
	//printf("x, R, B; %.5e\t%.5e\t%.5e\n", x, R, B); //matches

	    //----------------- B-field projection onto sky for this section --------------//
        // Now have introduced R dependence on the random blocks as well, works just like for the helical case, just * R_old/R_new on B_2 term
        //now can also include doppler depolarisation, changing the B field to EFFECTIVE B fields which give correct perp E component as if it had
        //been doppler depolarisedx
        //Now each block has its own theta_obs as we are looking inside the conical jet

        //Rotate B-fields along spherical cone surface slightly and also transversify depending on R/R0
        for (i=0; i<(n_blocks); i++){
            Bx[i] = B_0[i]*(cos(theta_r[i])+cos(M_PI/2+theta_phi[i])*cos(M_PI/2+theta_phi[i])*(1-cos(theta_r[i])))
                    + B_1[i]*cos(M_PI/2+theta_phi[i])*sin(M_PI/2+theta_phi[i])*(1-cos(theta_r[i]))
                    + B_2[i]*R0/R*sin(M_PI/2+theta_phi[i])*sin(theta_r[i]);
            By[i] = B_0[i]*sin(M_PI/2+theta_phi[i])*cos(M_PI/2+theta_phi[i])*(1-cos(theta_r[i]))
                    + B_1[i]*(cos(theta_r[i])+sin(M_PI/2+theta_phi[i])*sin(M_PI/2+theta_phi[i])*(1-cos(theta_r[i])))
                    - B_2[i]*R0/R*cos(M_PI/2+theta_phi[i])*sin(theta_r[i]);
            Bz[i] = -B_0[i]*sin(M_PI/2+theta_phi[i])*sin(theta_r[i])
                    + B_1[i]*cos(M_PI/2+theta_phi[i])*sin(theta_r[i])
                    + B_2[i]*R0/R*cos(theta_r[i]);
        //helical B-field equivalent has to be inside loop below as it depends on R
        }

        //Rotate B-fields along spherical cone surface slightly and also transversify depending on R/R0
        for (i=0; i<(n_blocks); i++){
            Bx_helical[i] = -R*sin(t_helix+phi_helix)*(cos(theta_r[i])+cos(M_PI/2+theta_phi[i])*cos(M_PI/2+theta_phi[i])*(1-cos(theta_r[i])))
                            + R*cos(t_helix+phi_helix)*cos(M_PI/2+theta_phi[i])*sin(M_PI/2+theta_phi[i])*(1-cos(theta_r[i]))
                            + c_helix*sin(M_PI/2+theta_phi[i])*sin(theta_r[i]);
            By_helical[i] = -R*sin(t_helix+phi_helix)*sin(M_PI/2+theta_phi[i])*cos(M_PI/2+theta_phi[i])*(1-cos(theta_r[i]))
                            + R*cos(t_helix+phi_helix)*(cos(theta_r[i])+sin(M_PI/2+theta_phi[i])*sin(M_PI/2+theta_phi[i])*(1-cos(theta_r[i])))
                            - c_helix*cos(M_PI/2+theta_phi[i])*sin(theta_r[i]);
            Bz_helical[i] = R*sin(t_helix+phi_helix)*sin(M_PI/2+theta_phi[i])*sin(theta_r[i])
                            + R*cos(t_helix+phi_helix)*cos(M_PI/2+theta_phi[i])*sin(theta_r[i])
                            + c_helix*cos(theta_r[i]);
        }

        /*REPLACE READING Bs into python with all STokes in C
        Now Find matrix of all B_effectives for all zones */
        if (!EVPA_rotation) {
            for (i=0; i<(n_blocks); i++) {
//                if (i==0){
//                    BX[i] = (Bz_helical[i]*sin(deg2rad(theta_obs)) + Bx_helical[i]*cos(deg2rad(theta_obs)));
//                    BY[i] = By_helical[i];
//                    BZ[i] = (Bz_helical[i]*cos(deg2rad(theta_obs))-Bx_helical[i]*sin(deg2rad(theta_obs)));
//                } else {
                    BX[i] = (Bz[i]*sin(deg2rad(theta_obs)) + Bx[i]*cos(deg2rad(theta_obs)));
                    BY[i] = By[i];
                    BZ[i] = (Bz[i]*cos(deg2rad(theta_obs))-Bx[i]*sin(deg2rad(theta_obs)));
//                }

            }
        } else {
            for (i=0; i<(n_blocks); i++) {
                if (i<19){
                    BX[i] = (Bz_helical[i]*sin(deg2rad(theta_obs)) + Bx_helical[i]*cos(deg2rad(theta_obs)));
                    BY[i] = By_helical[i];
                    BZ[i] = (Bz_helical[i]*cos(deg2rad(theta_obs))-Bx_helical[i]*sin(deg2rad(theta_obs)));
                } else {
                    BX[i] = (Bz[i]*sin(deg2rad(theta_obs)) + Bx[i]*cos(deg2rad(theta_obs)));
                    BY[i] = By[i];
                    BZ[i] = (Bz[i]*cos(deg2rad(theta_obs))-Bx[i]*sin(deg2rad(theta_obs)));
                }

            }
        }

        if (!DopDep){
            for (i=0; i<(n_blocks); i++) {
                for (l=0; l<n_blocks; l++){
                    Blength = sqrt(BX[i]*BX[i]+BY[i]*BY[i]+BZ[i]*BZ[i]); //normalization important to get correct z angle

                    B_effectives[l][i][0] = BX[i]/Blength;
                    B_effectives[l][i][1] = BY[i]/Blength;
                    B_effectives[l][i][2] = BZ[i]/Blength;

                }
            }
        } else {
            for (i=0; i<(n_blocks); i++) {
                for (l=0; l<n_blocks; l++){
                    DD_3Beffective(BX[i],BY[i],BZ[i],
                                cos(deg2rad(theta_obs))*sin(theta_r[l])*cos(theta_phi[l])+sin(deg2rad(theta_obs))*cos(theta_r[l]),
                                sin(theta_r[l])*sin(theta_phi[l]), -sin(deg2rad(theta_obs))*sin(theta_r[l])*cos(theta_phi[l])+cos(deg2rad(theta_obs))*cos(theta_r[l]),
                                gamma_bulk,B_effectives[l][i]);

                }
            }
        }

//        for (i=0; i<(n_blocks); i++) { //save B-fields for visualization
//            fprintf(Proj_Bfile,"\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n",theta_r[i],theta_phi[i],
//                    theta_tot[i],atan2(B_1[i],B_0[i]),atan2(By[i],Bx[i]),atan2(By_helical[i],Bx_helical[i]),
//                    atan2(BY[i],BX[i]),atan2(B_effectives[i][i][1],B_effectives[i][i][0]));
//        }

        //begin with Synchrotron emission
        for (i=0; i<array_size; i++)
        {
            P_single[i] = larmor(B, beta_e[i], gamma_e[i]); //gives the radiated power by one electron in each bin
            j[i] = j_per_hz(P_single[i], R, dN_dE[i], eps, f_c[i]);//*(R_next/R)*(R_next/R);//gives emissivity per unit frequency

            k[i] = k_new(j[i], eps, f_c[i]); //gives the corresponding opacity
            //Ps_per_m[i] = power_emitted(j[i], k[i], R)*M_PI*R*dfreqs[i]; //power per unit length assuming R roughly constant
            //Ps_per_m[i] = P_single[i]*(2*M_PI*pow(Me,2)*pow(C,2))*dN_dE[i]*dfreqs[i]/(C*3*gamma_e[i]*Qe*B);
            Ps_per_m_elec[i] = P_single[i]*dN_dE[i]*E_elecs[i]*Qe*dfreqs_pol[i] / (2*f_c[i]*C);
            //printf("Ps_per_m_elec %.5e\n", Ps_per_m_elec[i]);
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
                        P_perp[l] += (sqrt(3)*M_PI/(16*Me*pow(C,3)))*pow(Qe,3) * B * FG[m+1][1] * dN_dE[i]*dEe[i]*Qe * dfreqs_pol[l] * 4.4216E11 / (1E-7);     //power per unit frequency, have to *dfreq[j]
                        P_para[l] += (sqrt(3)*M_PI/(16*Me*pow(C,3)))*pow(Qe,3) * B * FG[m+1][2] * dN_dE[i]*dEe[i]*Qe * dfreqs_pol[l] * 4.4216E11 / (1E-7) ; //*dx/l_c !!! (this happens further down) and * other constants
                        //Ps_per_m_test[i] += (sqrt(3)*M_PI/(16*Me*pow(C,3)))*pow(Qe,3) * B * (FG[m+1][1]+FG[m+1][2]) * dN_dE[i]*dEe[i]*Qe * dfreqs_pol[l] * 4.4216E11 / (1E-7) ;


                    }
                    else if ((f_pol[l]/f_c[i]) <= FG[0][0] || (f_pol[l]/f_c[i]) >= FG[73][0])
                    {
                        P_perp[l] += 0.0;
                        P_para[l] += 0.0;
                        //Ps_per_m_test[i] += 0.0;
                    }
                }
            }
            Sync_losses[i] = Ps_per_m_elec[i];
            //printf("Ps_per_me_test %.5e\t%.5e\n", Ps_per_m_test[i],Ps_per_m_elec[i]);

        }

        for (l=0; l<array_size; l++){
            //P_para[l] = 0.2*P_perp[l];
            //make alpha dependent on frequency here
            if (l == 0){
               effective_alpha[l] = -(log10(dN_dE[l+1]) - log10(dN_dE[l]))/(log10(E_elecs[l+1]*Qe) - log10(E_elecs[l]*Qe));
            } else if (l > array_size - 2) {
                effective_alpha[l] = -(log10(dN_dE[l]) - log10(dN_dE[l-1]))/(log10(E_elecs[l]*Qe) - log10(E_elecs[l-1]*Qe));
            } else {
                effective_alpha[l] = -(log10(dN_dE[l+1]) - log10(dN_dE[l-1]))/(log10(E_elecs[l+1]*Qe) - log10(E_elecs[l-1]*Qe));
            }
            if (isnan(effective_alpha[l]) || isinf(effective_alpha[l])) {
                effective_alpha[l] = 8.41283e+01;
            }
            //printf("alphaef \t%.5e\n",effective_alpha[l]);
            for (g=0; g<n_blocks; g++) { //Sync Stokes Parameters
                q_theta = pow(sin(acos(B_effectives[g][g][2])),(alpha+1)/2); //factor to account for weaker emission when B_field pointed closer to line of sight
                S_Stokes[g][l][0] += q_theta*P_perp[l]/n_blocks, S_Stokes[g][l][1] += (q_theta*P_perp[l]/n_blocks)*cos(2*atan2(B_effectives[g][g][0],-B_effectives[g][g][1])), S_Stokes[g][l][2] += (q_theta*P_perp[l]/n_blocks)*sin(2*atan2(B_effectives[g][g][0],-B_effectives[g][g][1]));
                S_Stokes[g][l][0] += q_theta*P_para[l]/n_blocks, S_Stokes[g][l][1] += (q_theta*P_para[l]/n_blocks)*cos(2*atan2(B_effectives[g][g][1],B_effectives[g][g][0])), S_Stokes[g][l][2] += (q_theta*P_para[l]/n_blocks)*sin(2*atan2(B_effectives[g][g][1],B_effectives[g][g][0]));
                //printf("SSStokes123 %.5e\t%.5e\t%.5e\n", S_Stokes[g][l][0],S_Stokes[g][l][1],S_Stokes[g][l][2]);
            }

            Ps_per_m[l] = P_perp[l] + P_para[l];
            //printf("Ps_per_me_test, Ps_per_m %.5e\t%.5e\n", Ps_per_m_test[l],Ps_per_m[l]);
            //printf("Ps_permSYNC %.5e\n", Ps_per_m[l]);
            //photon energy density
            urad_array_perp[l] = 0.0; //initialise to avoid weird C bugs
            urad_array_para[l] = 0.0;
            //urad_array_perp[l] = P_perp[l]/(2.0*M_PI*C*R);
            //urad_array_perp[l] = P_perp[l]/(M_PI*R*R);
            urad_array_perp[l] = P_perp[l]/(2.0*M_PI*C*R);
            urad_array_para[l] = P_para[l]/(2.0*M_PI*C*R); //have these separate for compton polarisation calculation
            //u_rad += urad_array[l];
            //printf("Ne_e %.5e\n", Ne_e[l]);
            fprintf(photonpop, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", f_pol[l], (urad_array_perp[l]+urad_array_para[l]) /(f_pol[l]*H),(urad_array_perp[l]+urad_array_para[l]), f_em_max_IC[l], f_em_min_IC[l], f_pol_IC[l]);
        }
        //printf("step\n\n");
//        for (g=0; g<n_blocks; g++){
//            for (l=0; l<array_size; l++) {
//                printf("S_Stokes %.5e\n",0.5*atan(S_Stokes[g][l][2]/S_Stokes[g][l][1]));
//
//            }
//            printf("\n");
//
//        }



        //********************************************** ALP IC losses *******************************************//
        //assume Btheta fixed for now at 30deg, to avoid having to loop over every Bzone
        //double sig1=0.0;
        //double sig2=0.0;


//        for (g=0; g<n_blocks; g++){ //B-field blocks
//            for (h=0; h<n_blocks; h++){
//                Btheta = acos(B_effectives[h][g][2]); //compton polarization fraction very dependent on this for a single zone, 90deg gives highest
//                zeta = atan(B_effectives[h][g][1]/B_effectives[h][g][0]); //to rotate each Stokes to lab frame
//                //printf("Btheta [deg] \t%.5e", Btheta*180/M_PI);
//
//                if (h==g){
//                    zone_d = 1;
//                    phik_start = 0, phik_end = 10, cosk_start = 0, cosk_end = 10;
//                    //rotate to B_effective
//                    dfactor = 1;
//
//                } else {
//                    psi = - atan(BY[g]/BX[g]);
//                    zone_d = sqrt(pow(align[h][g][0]*cos(psi) - align[h][g][1]*sin(psi),2) + pow(align[h][g][0]*sin(psi) + align[h][g][1]*cos(psi),2) + pow(align[h][g][2],2));
//                    cosk_single = align[h][g][2]/zone_d;
//                    phik_single = atan2(align[h][g][0]*sin(psi) + align[h][g][1]*cos(psi), align[h][g][0]*cos(psi) - align[h][g][1]*sin(psi)); //this is between -pi and pi, but phik between 0 and 2pi
//
//
//                    cosk_start = findClosest(cosk, cosk_single, 10), cosk_end = cosk_start + 1; //findClosest might have a problem with 0.0, double check
//                    phik_start = findClosest(phik, phik_single, 10), phik_end = phik_start + 1;
//                    //printf("coskphik %.5e\t%d\t%.5e\t%d\n",cosk_single,cosk_start,phik_single,phik_start);
//                    //zone_d += 1;
//                    dfactor = 100 * asin(1/(2*zone_d))/M_PI;//this gives the 1/r dependence
//                    //printf("zoned %.5e\n", zone_d);
//                    //rotate to g's B_effective, but the B_effective given by the DOppler factor for that zone
//
//                }
//
//                for (n=0; n<array_size; n++){  //f_polIC (compton energy)
//                    for (l=0; l<array_size; l++){ //f_pol (sync energy)
//                        for (p=phik_start; p<phik_end; p++){ //phi
//                            for (m=cosk_start; m<cosk_end; m++){ //cosk
//                                F_min = sqrt(f_pol_IC[n]/(2*f_pol[l] * (1-cosk[m])));
////                                X = sqrt(f_pol_IC[n]/(2*f_pol[l]))*f_pol[l]*4.14E-15/Me_EV;
//                                v_k[0] = sqrt(1-cosk[m]*cosk[m])*cos(phik[p]), v_k[1] = sqrt(1-cosk[m]*cosk[m])*sin(phik[p]), v_k[2] = cosk[m]; //incoming photon direction vector
//                                e_k[0] = sqrt(1-cosk[m]*cosk[m])*sin(phik[p])*cos(Btheta);
//                                e_k[1] = sin(Btheta)*cosk[m] - cos(Btheta)*sqrt(1-cosk[m]*cosk[m])*cos(phik[p]); //incoming photon polarization vector (perp)
//                                e_k[2] = -sin(Btheta)*sqrt(1-cosk[m]*cosk[m])*sin(phik[p]);
//                                //this is just cross product of e_k and v_k above
//                                e_kpara[0] = pow(cosk[m],2)*sin(Btheta) - cos(Btheta)*cos(phik[p])*cosk[m]*sqrt(1-cosk[m]*cosk[m]) + (1-cosk[m]*cosk[m])*pow(sin(phik[p]),2)*sin(Btheta),
//                                e_kpara[1] = -sin(Btheta)*(1-cosk[m]*cosk[m])*sin(2*phik[p])/2 - cosk[m]*sqrt(1-cosk[m]*cosk[m])*sin(phik[p])*cos(Btheta), //incoming photon polarization vector (para)
//                                e_kpara[2] = (1-cosk[m]*cosk[m])*pow(sin(phik[p]),2)*cos(Btheta) - sqrt(1-cosk[m]*cosk[m])*cos(phik[p])*sin(Btheta)*cosk[m] + (1-cosk[m]*cosk[m])*pow(cos(phik[p]),2)*cos(Btheta);
//                                //e_k above is not normalised
//                                q_theta = pow(1-pow(cos(Btheta)*cosk[m]+sin(Btheta)*cos(phik[p])*sqrt(1-cosk[m]*cosk[m]),2),(effective_alpha[l]+1)/4);
//
//
////                                if (X > 1) {
////                                    KN = 0.4;
////                                } else if (X>10){
////                                    KN = 0.1;
////                                }
////                                else if (X>100){
////                                    KN = 0.02;
////                                }
////                                else if (X>1000){ // maybe too intense
////                                    KN = 0.0;
////                                }
////                                else {
////                                    KN = 1.0;
////                                }
//
//                                //printf("F_min,f_pol_IC[n],f_pol[l],cosk[m] %.5e\t%.5e\t%.5e\t%.5e\n", F_min,f_pol_IC[n],f_pol[l],cosk[m]);
//                                if (F_min*(Me_EV) > E_elecs[array_size-1]) {
//                                    P_perpIC[n] += 0.0;
//                                    P_paraIC[n] += 0.0;
//                                    continue;
//                                }
//                                else if (F_min*(Me_EV) > E_elecs[0]) {
//                                    for (o=0; o<array_size; o++){ //find Fmin location in E_elecs
//                                        if (F_min*(Me_EV) > E_elec_min[o] && F_min*(Me_EV) < E_elec_max[o]){
//                                        break;
//                                        }
//                                    }
//                                    for (i=o; i<array_size; i++){ //E_e
////                                        Pperpperp = 0.0;
////                                        Pperppara = 0.0;
////                                        Pparaperp = 0.0;
////                                        Pparapara = 0.0;
//
//                                        //3.95417e-103
//                                        //initial constant factor to make sure this matches with isotropic electron losses, original constant from paper is 1.2859E-91
//                                        Pperpperp = 3.95417e-103 /*8.26784E-118*/ * q_theta * dfreqs_polIC[n] * f_pol_IC[n]/f_pol[l] * F_min * /*dfreqs_pol[l]* */  (urad_array_perp[l] /(f_pol[l]*H)) * dcosk * dphik * dfactor
//                                        * ( Z_e(_e,e_k,v_k,1) * ( Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i])
//                                        + Sigma_1(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) ) + Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) )
//                                        *KN; //3/4*((1+X)/pow(X,3) * (2*X*(1+X)/(1+2*X) - log(1+2*X)) + 1/(2*X)*log(1+2*X) - (1+3*X)/pow(1+2*X,2)) ; //KN cross section correction for power cutoff (doesnt affect polarization)
//                                        //printf("Pperpperp %.5e\n", Pperpperp);
//                                        IC_Stokes[h][n][0] += (Pperpperp/(n_blocks)), IC_Stokes[h][n][1] += (Pperpperp/(n_blocks))*cos(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta))), IC_Stokes[h][n][2] += (Pperpperp/(n_blocks))*sin(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta)));
//
//                                        Pperppara = 3.95417e-103 /*8.26784E-118*/ * q_theta * dfreqs_polIC[n] * f_pol_IC[n]/f_pol[l] * F_min * /*dfreqs_pol[l]* */ (urad_array_para[l] /(f_pol[l]*H)) * dcosk * dphik * dfactor
//                                        * ( Z_e(_e,e_kpara,v_k,1) * ( Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i])
//                                        + Sigma_1(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) ) + Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) )
//                                        *KN; //3/4*((1+X)/pow(X,3) * (2*X*(1+X)/(1+2*X) - log(1+2*X)) + 1/(2*X)*log(1+2*X) - (1+3*X)/pow(1+2*X,2)) ; //KN cross section correction for power cutoff (doesnt affect polarization)
//                                        //printf("Pperppara %.5e\n", Pperppara);
//                                        IC_Stokes[h][n][0] += (Pperppara/(n_blocks)), IC_Stokes[h][n][1] += (Pperppara/(n_blocks))*cos(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta))), IC_Stokes[h][n][2] += (Pperppara/(n_blocks))*sin(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta)));
//                                        //printf("ICStokes123 %.5e\t%.5e\t%.5e\n", IC_Stokes[h][n][0],IC_Stokes[h][n][1],IC_Stokes[h][n][2]);
//
//                                        Pparaperp = 3.95417e-103 /*8.26784E-118*/ * q_theta * dfreqs_polIC[n] * f_pol_IC[n]/f_pol[l] * F_min * /*dfreqs_pol[l]* */ (urad_array_perp[l] /(f_pol[l]*H)) * dcosk * dphik * dfactor
//                                        * ( Z_e(_e,e_k,v_k,0) * ( Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i])
//                                        + Sigma_1(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) ) + Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) )
//                                        *KN; //3/4*((1+X)/pow(X,3) * (2*X*(1+X)/(1+2*X) - log(1+2*X)) + 1/(2*X)*log(1+2*X) - (1+3*X)/pow(1+2*X,2)) ; //KN cross section correction for power cutoff (doesnt affect polarization)
//                                        //printf("Pparaperp %.5e\n", Pperpperp);
//                                        IC_Stokes[h][n][0] += (Pparaperp/(n_blocks)), IC_Stokes[h][n][1] += (Pparaperp/(n_blocks))*cos(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta))), IC_Stokes[h][n][2] += (Pparaperp/(n_blocks))*sin(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta)));
//
//                                        //printf("Pparaperp1 %.5e\n", Pparaperp);
//
//                                        Pparapara = 3.95417e-103 /*8.26784E-118*/ * q_theta * dfreqs_polIC[n] * f_pol_IC[n]/f_pol[l] * F_min * /*dfreqs_pol[l]* */ (urad_array_para[l] /(f_pol[l]*H)) * dcosk * dphik * dfactor
//                                        * ( Z_e(_e,e_kpara,v_k,0) * ( Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i])
//                                        + Sigma_1(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) ) + Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) )
//                                        *KN; //3/4*((1+X)/pow(X,3) * (2*X*(1+X)/(1+2*X) - log(1+2*X)) + 1/(2*X)*log(1+2*X) - (1+3*X)/pow(1+2*X,2)) ; //KN cross section correction for power cutoff (doesnt affect polarization)
//                                        //printf("Pparapara %.5e\n", Pparapara);
//                                        IC_Stokes[h][n][0] += (Pparapara/(n_blocks)), IC_Stokes[h][n][1] += (Pparapara/(n_blocks))*cos(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta))), IC_Stokes[h][n][2] += (Pparapara/(n_blocks))*sin(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta)));
//
//                                        P_perpIC[n] += Pperpperp + Pperppara;
//                                        P_paraIC[n] += Pparaperp + Pparapara;
//
//                                        //Ps_per_m_test[i] += Pperpperp + Pperppara + Pparaperp + Pparapara;
//                                    }
//                                    //for (i=o; i<array_size; i++){
//                                    //    sig1 += Sigma_1(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]);
//                                    //    sig2 += Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]);
//                                    //}
//                                    //printf("PiBCS %.5e\t%.5e\t%.5d\t%.5e\t%.5e\n", (sig1+sig2)/(sig1+3*sig2), F_min,1,sig1,sig2);
//                                    //sig1 = 0.0;
//                                    //sig2 = 0.0;
//
//                                }
//                                else {
//                                    for (i=0; i<array_size; i++){ //E_e
////                                        Pperpperp = 0.0;
////                                        Pperppara = 0.0;
////                                        Pparaperp = 0.0;
////                                        Pparapara = 0.0;
//
//                                        Pperpperp = 3.95417e-103 /*8.26784E-118*/ * q_theta * dfreqs_polIC[n] * f_pol_IC[n]/f_pol[l] * F_min * /*dfreqs_pol[l]* */  (urad_array_perp[l] /(f_pol[l]*H)) * dcosk * dphik * dfactor
//                                        * ( Z_e(_e,e_k,v_k,1) * ( Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i])
//                                        + Sigma_1(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) ) + Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) )
//                                        *KN; //3/4*((1+X)/pow(X,3) * (2*X*(1+X)/(1+2*X) - log(1+2*X)) + 1/(2*X)*log(1+2*X) - (1+3*X)/pow(1+2*X,2)) ; //KN cross section correction for power cutoff (doesnt affect polarization)
//                                        //printf("Pperpperp1 %.5e\n", Pperpperp);
//                                        IC_Stokes[h][n][0] += (Pperpperp/(n_blocks)), IC_Stokes[h][n][1] += (Pperpperp/(n_blocks))*cos(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta))), IC_Stokes[h][n][2] += (Pperpperp/(n_blocks))*sin(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta)));
//
//                                        Pperppara = 3.95417e-103 /*8.26784E-118*/ * q_theta * dfreqs_polIC[n] * f_pol_IC[n]/f_pol[l] * F_min * /*dfreqs_pol[l]* */ (urad_array_para[l] /(f_pol[l]*H)) * dcosk * dphik * dfactor
//                                        * ( Z_e(_e,e_kpara,v_k,1) * ( Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i])
//                                        + Sigma_1(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) ) + Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) )
//                                        *KN; //3/4*((1+X)/pow(X,3) * (2*X*(1+X)/(1+2*X) - log(1+2*X)) + 1/(2*X)*log(1+2*X) - (1+3*X)/pow(1+2*X,2)) ; //KN cross section correction for power cutoff (doesnt affect polarization)
//                                        //printf("Pperppara1 %.5e\n", Pperppara);
//                                        IC_Stokes[h][n][0] += (Pperppara/(n_blocks)), IC_Stokes[h][n][1] += (Pperppara/(n_blocks))*cos(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta))), IC_Stokes[h][n][2] += (Pperppara/(n_blocks))*sin(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta)));
//                                        //printf("ICStokes123 %.5e\t%.5e\t%.5e\n", IC_Stokes[h][n][0],IC_Stokes[h][n][1],IC_Stokes[h][n][2]);
//
//                                        Pparaperp = 3.95417e-103 /*8.26784E-118*/ * q_theta * dfreqs_polIC[n] * f_pol_IC[n]/f_pol[l] * F_min * /*dfreqs_pol[l]* */ (urad_array_perp[l] /(f_pol[l]*H)) * dcosk * dphik * dfactor
//                                        * ( Z_e(_e,e_k,v_k,0) * ( Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i])
//                                        + Sigma_1(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) ) + Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) )
//                                        *KN; //3/4*((1+X)/pow(X,3) * (2*X*(1+X)/(1+2*X) - log(1+2*X)) + 1/(2*X)*log(1+2*X) - (1+3*X)/pow(1+2*X,2)) ; //KN cross section correction for power cutoff (doesnt affect polarization)
//                                        //printf("Pparaperp1 %.5e\n", Pperpperp);
//                                        IC_Stokes[h][n][0] += (Pparaperp/(n_blocks)), IC_Stokes[h][n][1] += (Pparaperp/(n_blocks))*cos(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta))), IC_Stokes[h][n][2] += (Pparaperp/(n_blocks))*sin(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta)));
//
//                                        //printf("Pparaperp1 %.5e\n", Pparaperp);
//
//                                        Pparapara = 3.95417e-103 /*8.26784E-118*/ * q_theta * dfreqs_polIC[n] * f_pol_IC[n]/f_pol[l] * F_min * /*dfreqs_pol[l]* */ (urad_array_para[l] /(f_pol[l]*H)) * dcosk * dphik * dfactor
//                                        * ( Z_e(_e,e_kpara,v_k,0) * ( Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i])
//                                        + Sigma_1(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) ) + Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]) )
//                                        *KN; //3/4*((1+X)/pow(X,3) * (2*X*(1+X)/(1+2*X) - log(1+2*X)) + 1/(2*X)*log(1+2*X) - (1+3*X)/pow(1+2*X,2)) ; //KN cross section correction for power cutoff (doesnt affect polarization)
//                                        //printf("Pparapara1 %.5e\n", Pparapara);
//                                        IC_Stokes[h][n][0] += (Pparapara/(n_blocks)), IC_Stokes[h][n][1] += (Pparapara/(n_blocks))*cos(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta))), IC_Stokes[h][n][2] += (Pparapara/(n_blocks))*sin(2*atan2(_e[1]*cos(zeta)+_e[0]*sin(zeta),_e[0]*cos(zeta)-_e[1]*sin(zeta)));
//
//                                        //printf("Pparapara2 %.5e\n", Pparapara);
//
//                                        P_perpIC[n] += Pperpperp + Pperppara;
//                                        P_paraIC[n] += Pparaperp + Pparapara;
//
//                                        //Ps_per_m_test[i] += Pperpperp + Pperppara + Pparaperp + Pparapara;
//                                    }
//                                    //for (i=0; i<array_size; i++){
//                                    //    sig1 += Sigma_1(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]);
//                                    //    sig2 += Sigma_2(E_elecs[i],Me_EV*Qe*F_min,dN_dE[i],dEe[i]);
//                                    //}
//                                    //printf("PiBCS %.5e\t%.5e\t%.5d\t%.5e\t%.5e\n", (sig1+sig2)/(sig1+3*sig2),F_min,2,sig1,sig2);
//                                    //sig1 = 0.0;
//                                    //sig2 = 0.0;
//                                }
//                            }
//                        }
//                    }
//                }
//
//            }
////            for (n=0; n<array_size; n++){
////                fprintf(PperpfileIC, "\t%.5e", P_perpIC[n]);//*pow(doppler_factor, 4.0)); //polarisation power now doppler boosted in python file
////                fprintf(PparafileIC, "\t%.5e",P_paraIC[n]);//*pow(doppler_factor, 4.0));
////                if (n==array_size-1) {
////                    fprintf(PperpfileIC, "\n");
////                    fprintf(PparafileIC, "\n");
////                }
////            }
//        }


        //printf("ICPI %.5e\t%.5e\t%.5e\n", IC_Pi[30],IC_PA[10],S_Pi[45]);

        for (i=0; i<array_size; i++){
            u_rad += (urad_array_perp[i]+urad_array_para[i]);
            //printf("Ps_permtest %.5e\n", Ps_per_m_test[i]);
            //printf("u_rad %.5e\n", u_rad);
        }
        
        for (i=0; i<array_size; i++) //1.79 //0.85
        {
            Ps_per_m_IC[i] = 0.0;
            Ps_per_mIC_elec[i] = 0.0;
            Ps_per_m_IC[i] = P_perpIC[i] + P_paraIC[i];
            //printf("ps %.5e\n", Ps_per_m_IC[i]);
            Ps_per_mIC_elec[i] = (4/3)*6.6524E-29*gamma_e[i]*gamma_e[i]*beta_e[i]*beta_e[i]*dN_dE[i]*dEe[i]*Qe*u_rad; //test for Ps_per_mIC_elec == this
            IC_losses[i] = Ps_per_mIC_elec[i];
            //printf("ps_elecIC \t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n", gamma_e[i],beta_e[i],dN_dE[i],dEe[i],u_rad,Qe);

        }
        
        
        //**********************************************************************************************//
        
        
        
        //some code to estimate what dx should be
        dx_set = findminelementNO0(Ne_e, array_size);//gets the lowest non_zero element. Higher energy electrons radiate more rapidly


        dx_R = 0.05*(R0+x*tan(deg2rad(theta_open_p)))/tan(deg2rad(theta_open_p)); //ensures Rnew <= 1.05 Rold
        dx_P = (Ne_e[dx_set]*(E_elecs[dx_set]-E_elecs[dx_set-1])*Qe)/(Ps_per_m_elec[dx_set] + Ps_per_mIC_elec[dx_set]); //based on radiative losses
        //printf(" dx_P  %.5e\t%d \n", Ps_per_mIC_elec[dx_set],dx_set);
        //printf("dx_R, dx_P \t%.5e\t%.5e \n", dx_P, (Ne_e[dx_set]*(E_elecs[dx_set]-E_elecs[dx_set-1])*Qe)/(Ps_per_m[dx_set]));
  
        if (dx_R < dx_P)
        {
            dx = dx_R;
        }

        if (dx_P < dx_R)
        {
            dx = dx_P;
        }

        //printf("dx_R dx_P dx %.5e\t%.5e\t%.5e \n", dx_R, dx_P, dx);

        for (h=0; h<n_blocks; h++){ //summing the Stokes parameters of all the blocks and boosting
            zone_doppler = doppler(beta_bulk, theta_tot[h]);
            //difference between bin boundaries in logspace:
            logdif = log10(f_pol_IC[1])-log10(f_pol_IC[0]);
            binshiftIC = (int)round(log10(zone_doppler/doppler_factor)/logdif); //the number of bins one moves across (relative to shift of bins which already takes place)
            logdif = log10(f_pol[1])-log10(f_pol[0]);
            binshiftS = (int)round(log10(zone_doppler/doppler_factor)/logdif);
            printf("binshifts %d\t%d\n", binshiftIC,binshiftS);
            printf("dopplers %.5e\t%.5e\t%.5e\n", zone_doppler,doppler_factor,logdif);
            //binshiftS = 0;
            //binshiftIC = 0;
            //think about how binsizes change

            for (n=0; n<array_size; n++){
                //printf("ICStokes123 %.5e\t%.5e\t%.5e\n", IC_Stokes[h][n][0],IC_Stokes[h][n][1],IC_Stokes[h][n][2]);
                if (n-binshiftIC >= 0 && n-binshiftIC < array_size){
                    IC_StokesTotal[n][0] += IC_Stokes[h][n-binshiftIC][0] * pow(zone_doppler,4) * dx;//(dfreqs_polIC[n-binshiftIC]/dfreqs_polIC[n]) * dx;
                    IC_StokesTotal[n][1] += IC_Stokes[h][n-binshiftIC][1] * pow(zone_doppler,4) * dx;//(dfreqs_polIC[n-binshiftIC]/dfreqs_polIC[n]) * dx;
                    IC_StokesTotal[n][2] += IC_Stokes[h][n-binshiftIC][2] * pow(zone_doppler,4) * dx;//(dfreqs_polIC[n-binshiftIC]/dfreqs_polIC[n]) * dx;
                }

                if (n-binshiftS >= 0 && n-binshiftS < array_size){
                    //printf("SSStokes123 %.5e\t%.5e\t%.5e\n", S_Stokes[h][n][0],S_Stokes[h][n][1],S_Stokes[h][n][2]);
                    S_StokesTotal[n][0] += S_Stokes[h][n-binshiftS][0] * pow(zone_doppler,4) * dx; //(dfreqs_pol[n-binshiftS]/dfreqs_pol[n]) * dx;
                    S_StokesTotal[n][1] += S_Stokes[h][n-binshiftS][1] * pow(zone_doppler,4) * dx; //(dfreqs_pol[n-binshiftS]/dfreqs_pol[n]) * dx;
                    S_StokesTotal[n][2] += S_Stokes[h][n-binshiftS][2] * pow(zone_doppler,4) * dx; //(dfreqs_pol[n-binshiftS]/dfreqs_pol[n]) * dx;
                }
                //printf("totALS \t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n", S_Stokes[h][n][0],S_Stokes[h][n][1],S_Stokes[h][n][2],IC_Stokes[h][n][0],IC_Stokes[h][n][1],IC_Stokes[h][n][2]);
                //printf("doppler \t%.5e\n", pow(zone_doppler,4));
            }
        }
        //printf("ICStokes123 %.5e\t%.5e\t%.5e\n", IC_StokesTotal[5][0],IC_StokesTotal[15][1],IC_StokesTotal[35][2]);

//        if (!nSteps){
//            for (n=0; n<array_size; n++){
//                IC_Pi[n] = sqrt(IC_StokesTotal[n][1]*IC_StokesTotal[n][1] + IC_StokesTotal[n][2]*IC_StokesTotal[n][2]) / IC_StokesTotal[n][0];
//                S_Pi[n] = sqrt(S_StokesTotal[n][1]*S_StokesTotal[n][1] + S_StokesTotal[n][2]*S_StokesTotal[n][2]) / S_StokesTotal[n][0];
//
//                IC_PA[n] = 0.5*atan(IC_StokesTotal[n][2]/IC_StokesTotal[n][1]);
//                S_PA[n] = 0.5*atan(S_StokesTotal[n][2]/S_StokesTotal[n][1]);
//
//                IC_P[n] = IC_StokesTotal[n][0];
//                S_P[n] = S_StokesTotal[n][0];
//
//                fprintf(TESTFIL, "\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n",S_Pi[n],S_PA[n],S_P[n],IC_Pi[n],IC_PA[n],IC_P[n]);
//            }
//        }
        
        //change values for next loop
        x += dx;
        R_prev = R;
        R = R_new(R0, deg2rad(theta_open_p), x);

        //correct synchrotron for full length + IC
        for (i=0; i<array_size; i++)
        {
            Ps_per_m[i] *= dx;
            Ps_per_m_IC[i] *= dx;
            Sync_losses[i] *= dx;
            IC_losses[i] *= dx;              //ALP
            P_perp[i] *= dx;  //polarisation powers
            P_para[i] *= dx;
            P_perpIC[i] *= dx;
            P_paraIC[i] *= dx;
	        fprintf(jfile, "\t%.5e", j[i]);
            fprintf(freqfile, "\t%.5e", f_c[i]); //these are Doppler boosted later
            fprintf(opfile, "\t%.5e", k[i]);//*(R/R_prev)*(R/R_prev));
            fprintf(powfile, "\t%.5e", Ps_per_m[i]*pow(doppler_factor, 4.0));//*pow((doppler(beta_e[i], deg2rad(theta_obs))), 4.0));//*Ne_e[i]);
            fprintf(ICfile, "\t%.5e", Ps_per_m_IC[i]*pow(doppler_factor, 4.0)); //ALP
            fprintf(Pperpfile, "\t%.5e", P_perp[i]);//*pow(doppler_factor, 4.0)); //polarisation power now doppler boosted in python file
            fprintf(Pparafile, "\t%.5e",P_para[i]);//*pow(doppler_factor, 4.0));
            fprintf(PperpfileIC, "\t%.5e", P_perpIC[i]);//*pow(doppler_factor, 4.0));
            fprintf(PparafileIC, "\t%.5e",P_paraIC[i]);//*pow(doppler_factor, 4.0));
            //fprintf(TESTFIL, "\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n",S_Pi[i],S_PA[i],S_P[i],IC_Pi[i],IC_PA[i],IC_P[i]);
            //fprintf(Proj_Bfile,"\t%.5e",Proj_theta_B);
            //fprintf(betadata, "\t%.5e", beta_e[i]);

            if (i==array_size-1) //puts the outputs for the next jet section on a new line in files
            {
                fprintf(freqfile, "\n");
                fprintf(opfile, "\n");
        		fprintf(jfile, "\n");
                fprintf(powfile, "\n");
                fprintf(ICfile, "\n");
                fprintf(Pperpfile, "\n");
                fprintf(Pparafile, "\n");
                fprintf(PperpfileIC, "\n");
                fprintf(PparafileIC, "\n");

                //fprintf(Proj_Bfile, "\n");
            }
        }
        //same thing as above but now just to save Proj_Bs as this has a different size (nblocks)
//        for(i=0; i<n_blocks; i++){
//            //save Projected B field angle for each jet section and block here
//            fprintf(Proj_Bfile,"\t%.5e",Proj_theta_B[i]);
//            if (i==n_blocks-1){
//                fprintf(Proj_Bfile,"\n");
//            }
//        }



        //compute the synchrotron losses + IC LOSSES
        for (i=0; i<array_size; i++)
        {
            if (i==0)
            {
                //Ne_losses[i]=0.0; //cannot drop to lower bins
                Ne_losses[i]=(Sync_losses[i] + IC_losses[i])*(C/C)/((E_elecs[i]-Me_EV)*Qe); //can't radiate anymore
            }
            else
            {
	      //Ne_losses[i] = (Sync_losses[i])*(C/C)/((E_elecs[i]-E_elecs[i-1])*Qe);
	      Ne_losses[i] = (Sync_losses[i] + IC_losses[i])*(C/C)/((E_elecs[i]-E_elecs[i-1])*Qe);
            }


            if (i==array_size-1)
            {
                Ne_gains[i] = 0.0;
            }
            else
            {
	      //Ne_gains[i] = (Sync_losses[i+1])*(C/C)/((E_elecs[i+1]-E_elecs[i])*Qe);
	      Ne_gains[i] = (Sync_losses[i+1] + IC_losses[i+1])*(C/C)/((E_elecs[i+1]-E_elecs[i])*Qe);
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
            dN_dE[i] = electron_PL(A_elecs[i], alpha, E_elecs[i]*Qe, E_max*Qe);
            P_para[i] = 0.0;//Reset
            P_perp[i] = 0.0;
            P_paraIC[i] = 0.0;
            P_perpIC[i] = 0.0;
            Ps_per_mIC_elec[i] = 0.0;
            Ps_per_m_elec[i] = 0.0;

            fprintf(elecpopfull, "%.6e\t%.6e\t%.6e\n", E_elecs[i], Ne_e[i], dN_dE[i]);

        }
        memset(S_Stokes, 0, sizeof(S_Stokes[0][0][0])* n_blocks * array_size * 3); //reset these to 0 for next section
        memset(IC_Stokes, 0, sizeof(IC_Stokes[0][0][0])* n_blocks * array_size * 3);



        nSteps+=1; //allows the number of jet sections to be determined
        //printf("nStep, dx, x, B, R %.5d\t%.5e\t%.5e\t%.5e\t%.5e\n", nSteps, dx, x, B, R);


        //save data needed at every jet section to solve line of sight opacity
        fprintf(basicdata, "\t%.5e\t%.5e\t%.5e\t%.5e \n", dx, x, B, R);

        //break;

    }

    int n_xray, n_opt, n_radio; //need dfactor to account for frequency boost
    n_xray = findClosest(f_pol,9.67E17/doppler_factor , array_size); //X-ray band 4keV
    n_opt = findClosest(f_pol,6E14/doppler_factor , array_size); //Optical band 2.5 eV
    n_radio = findClosest(f_pol,1.2E11/doppler_factor , array_size); //Radio band 5e-4 eV
    //printf("x,o,r %d\t%d\t%d\n",n_xray,n_opt,n_radio);


    for (n=0; n<array_size; n++){
            //IC_Pi[n] = sqrt(IC_StokesTotal[n][1]*IC_StokesTotal[n][1] + IC_StokesTotal[n][2]*IC_StokesTotal[n][2]) / IC_StokesTotal[n][0];
            S_Pi[n] = sqrt(S_StokesTotal[n][1]*S_StokesTotal[n][1] + S_StokesTotal[n][2]*S_StokesTotal[n][2]) / S_StokesTotal[n][0];

            //IC_PA[n] = 0.5*atan2(IC_StokesTotal[n][2],IC_StokesTotal[n][1]);
            S_PA[n] = 0.5*atan2(S_StokesTotal[n][2],S_StokesTotal[n][1]);

            //IC_P[n] = IC_StokesTotal[n][0];
            S_P[n] = S_StokesTotal[n][0];

            fprintf(TESTFIL2, "\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n",S_Pi[n],S_PA[n],
                   S_P[n],IC_Pi[n],IC_PA[n],IC_P[n]);

    }
//    fprintf(TESTFIL, "\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n",S_Pi[n_radio],S_Pi[n_opt],
//            S_Pi[n_xray],S_PA[n_radio],S_PA[n_opt],S_PA[n_xray]);



    printf("The jet was divided into %d %s \n", nSteps, "section(s).");

    //store electron information
    for (i=0; i<array_size; i++)

    {
        fprintf(elecpop, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e \n", E_elecs[i], Ne_orig[i], Ne_intermed[i], Ne_e[i], dN_dE_orig[i]); //file elecpop.txt

    }

    return 0;
}


