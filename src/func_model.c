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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>

#define Me_EV 0.511E6 //electron rest energy
#define C 3.0E8       //speed of light
#define Qe 1.6E-19    //electron charge
#define Me 9.11E-31   //electron mass kg
#define H 6.63E-34    //Plancks constant

#define STR_PRINT(x) #x
#ifdef DEBUG
#define PRINT(message, parameter)   \
    do                              \
    {                               \
        printf(message, parameter); \
    } while (0)
#else
#define PRINT(message, parameter)
#endif

const size_t ARRAY_SIZE = 50;        // sets the number of synchrotron & SSC bins
const double ARRAY_SIZE_D = 50.0; // use to set log ratios, should be the same as above line but with .0
size_t dx_set; // use to define smallest non zero population
// #ifdef DEBUG
// int SSC = 0; // SSC off for debugging
// #else
size_t SSC = 1; //SSC on by default
size_t nLT_sections = 30;
// #endif


int jetmodel(double *argv, int *argblocks, double *IC_StokesTotal, double *S_StokesTotal, double *f_pol_IC, double *f_pol) //argc is integer number of arguments passed, argv[0] is program name, argv[1..n] are arguments passed in string format
{
    int i, l, m, n, o, p, g, h;   //some looping parameters
    double c, d, q;

    // Jet Parameters
    double W_j = argv[0];             //1.3E37; // W jet power in lab frame. Should be OBSERVED POWER
    double L_jet = 5E20;              // length in m in the fluid frame
    double E_min = argv[8];           //5.11E6; // Minimum electron energy 
    double E_max = argv[1];           //1.7E10; // Energy of the ECO in eV 
    double alpha = argv[2];           //1.85; // PL index of electrons
    double theta_open_p = argv[3];    //40.0; // opening angle of the jet in the fluid frame 
    double gamma_bulk = argv[4];      //16; // bulk Lorentz factor of jet material
    double B = argv[5], B0 = argv[5]; //1E-4; // B-field at jet base
    double R0 = 0.0, R = 0.0;         // Radius of the jet at the base 3.32 works fairly well 
    double B_prev = 0.0;              // changing parameters of the jet-initialise. R prev corrects for increasing jet volume 
    double theta_obs = argv[6];       // observers angle to jet axis in rad 
    double A_eq = argv[7];            //1.0
    const size_t N_BLOCKS = argblocks[0];      // for the TEMZ model, can have 1,7,19,37,61,91,127 blocks, (rings 0,1,2,3,4,5,6)
    const size_t N_RINGS = argblocks[1];       // up to 6 rings possible atm, must choose number of rings corresponding to number of zones
    printf("Jet Parameters: %.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n", W_j, E_max, alpha, theta_open_p, gamma_bulk, B0, theta_obs, A_eq, E_min);

    // Define radius at jet base
    R0 = R_0(W_j * (3.0 / 4.0), A_eq, gamma_bulk, B0); //assumes equipartition fraction 1.0 usually
    R = R_0(W_j * (3.0 / 4.0), A_eq, gamma_bulk, B0);
    double R_prev = R0; //initialize
    printf("Radius at jet base: %.5e \n", R);

    // Read in any required data (only Bessel Functions for now)
    FILE *FGfile;
    FGfile = fopen("src/FG.txt", "r"); //Bessel Functions for Synchrotron

    double FG[75][3]; //Bessel Functions F, G and FG
    for (m = 0; m < 75; m++)
    {
        fscanf(FGfile, "%lf\t%lf\t%lf\n", &q, &c, &d);
        FG[m][0] = q;
        FG[m][1] = c;
        FG[m][2] = d;
    }

    //define some useful parameters to gauge progress
    size_t nSteps = 0; //how many sections the jet has broken down into
    double x = 0;   //progress along the jet
    double dx = 0;  //incremenet of step length
    double eps;     //stores the eqn 2.18b
    double beta_bulk = lortovel(gamma_bulk);
    double doppler_factor = doppler(beta_bulk, deg2rad(theta_obs));

    //define some parameters for determining the dx of each jet section
    double dx_R; //dx where R_new = 1.05*R_old
    double dx_P; //dx to ensure smallest bin population not depleted

    //some parameters used to define the population of electrons
    double A_elecs[ARRAY_SIZE];     //stores PL coefficient
    double Ne_e[ARRAY_SIZE];        // no of electrons in each bin
    double Ne_orig[ARRAY_SIZE];     // at jet base: allows comparison later
    double dN_dE[ARRAY_SIZE];       // no of electrons per J in each bin
    double dN_dE_orig[ARRAY_SIZE];
    double dEe[ARRAY_SIZE];       //size of bin widths in eV
    double Ne_losses[ARRAY_SIZE]; //losses per section
    double Ne_gains[ARRAY_SIZE];  //gains per section
    double A_orig;                // PL prefactor initially the same for all electrons
    double Ue_jetbase = 0.0;      // energy in particles at conical region base
    double UB_jetbase = 0.0;      //energy in particles at the jet base

    // parameters used to initialise electron population
    double Ee_max = 4.0E12;                                   //Me_EV*3.0E4;//4.0E12; //sets the maximum electron energy in eV
    double ratio = pow((Ee_max / E_min), (1 / ARRAY_SIZE_D)); //logarithmic ratio of electron bins: USED FOR BOUNDS
    double E_bounds[ARRAY_SIZE + 1], indices[ARRAY_SIZE + 1];
    arange(indices, ARRAY_SIZE + 1); //as np.arange

    for (i = 0; i < ARRAY_SIZE + 1; i++)
    {
        E_bounds[i] = E_min * (pow(ratio, indices[i]));
    }

    //now sort the bounds into max, min and E_elecs
    double E_elecs[ARRAY_SIZE];                            // stores the energy in eV of each electron
    double E_elec_min[ARRAY_SIZE], E_elec_max[ARRAY_SIZE]; //these map onto critical frequencies, should be in eV

    for (i = 0; i < ARRAY_SIZE; i++)
    {
        E_elec_min[i] = E_bounds[i];
        E_elec_max[i] = E_bounds[i + 1];
        E_elecs[i] = pow(E_elec_min[i] * E_elec_max[i], 0.5);
    }

    //parameters based on individual electron bins
    double gamma_e[ARRAY_SIZE]; // store Lorentz factors
    double beta_e[ARRAY_SIZE];  // v/c

    //parameters used in determining synchrotron emission
    double f_c[ARRAY_SIZE];                                                                              //store the critical frequencies
    double dfreqs[ARRAY_SIZE];                                                                           // store the bin widths in critical frequency
    double f_em_min[ARRAY_SIZE], f_em_max[ARRAY_SIZE], f_em_min_IC[ARRAY_SIZE], f_em_max_IC[ARRAY_SIZE]; // store bin boundaries for critical frequencies
    double j[ARRAY_SIZE];                                                                                // emissivity per Hz per unit volume
    double k[ARRAY_SIZE];                                                                                // opacity, units of m^-1
    double Ps_per_m[ARRAY_SIZE];                                                                         // Synchrotron power per unit m
    double Ps_per_m_elec[ARRAY_SIZE];                                                                    //power lost by each electron bin per m
    double Sync_losses[ARRAY_SIZE];     // store total electron energy losses
    double P_single[ARRAY_SIZE];        // power emitted by a single electron of energy E
    double Ps_per_mIC_elec[ARRAY_SIZE]; //energy lost by each electron bin per m of IC emission
    memset(Ps_per_mIC_elec, 0.0, ARRAY_SIZE * sizeof(Ps_per_mIC_elec[0]));
    double IC_losses[ARRAY_SIZE]; //store electron IC energy losses

    //******************************** Inititalise seed photon energy density variables *************************************//

    //want to allocate these as dynamic arrays so can change their size with realloc
    double *R_array, *dx_array, **Pperp_array, **Ppara_array;
    double *ptri, **ptrii;
    double counter = 0.0;
    double dxsum = 0;
    double Area0;
    double dist = 0;

    double urad_array_perp[N_BLOCKS][ARRAY_SIZE]; //photon energy density
    memset(urad_array_perp, 0.0, N_BLOCKS * ARRAY_SIZE * sizeof(urad_array_perp[0][0]));
    double urad_array_para[N_BLOCKS][ARRAY_SIZE];
    memset(urad_array_para, 0.0, N_BLOCKS * ARRAY_SIZE * sizeof(urad_array_para[0][0]));
    double dfactor_perp[N_BLOCKS][N_BLOCKS][ARRAY_SIZE]; //photon energy density
    memset(dfactor_perp, 0.0, N_BLOCKS * N_BLOCKS * ARRAY_SIZE * sizeof(dfactor_perp[0][0][0]));
    double dfactor_para[N_BLOCKS][N_BLOCKS][ARRAY_SIZE];
    memset(dfactor_para, 0.0, N_BLOCKS * N_BLOCKS * ARRAY_SIZE * sizeof(dfactor_para[0][0][0]));
    
    int buffer_subset[N_BLOCKS]; //different blocks have different buffer_sizes depending on their position in jet cross section
    for (i = 0; i < N_BLOCKS; i++)
    { //cant allocate non zero with memset
        buffer_subset[i] = 1;
    }
    int buffer_size = 1;

    //****************Initialise Polarisation variables *****************************************************************************//
    //first set frequency bins for polarised powers to drop into:
    double P_perp[ARRAY_SIZE];
    memset(P_perp, 0.0, ARRAY_SIZE * sizeof(P_perp[0]));
    double P_para[ARRAY_SIZE];
    memset(P_para, 0.0, ARRAY_SIZE * sizeof(P_para[0])); //initialise to make sure junk values arent inside upon first +=
    double dfreqs_pol[ARRAY_SIZE]; //again fixed size
    double dfreqs_polIC[ARRAY_SIZE];

    //Fill electron arrays with initial parameters
    for (i = 0; i < ARRAY_SIZE; i++)
    {
        E_elec_min[i] = E_elecs[i] * (1.0 / pow(ratio, 0.5));                    //lower electron eV bounds
        E_elec_max[i] = E_elecs[i] * pow(ratio, 0.5);                            //upper electron eV bounds
        dEe[i] = E_elec_max[i] - E_elec_min[i];                                  // sets bin widths in eV
        A_elecs[i] = A_PL(alpha, W_j, gamma_bulk, E_min * Qe, E_max * Qe, A_eq); //1 assumes equipartition. DIVIDED W_j by gamma^2 for jet frame
        dN_dE[i] = electron_PL(A_elecs[i], alpha, E_elecs[i] * Qe, E_max * Qe);
        dN_dE_orig[i] = electron_PL(A_elecs[i], alpha, E_elecs[i] * Qe, E_max * Qe);
        Ne_e[i] = electron_PL(A_elecs[i], alpha, E_elecs[i] * Qe, E_max * Qe) * dEe[i] * Qe;
        Ne_orig[i] = electron_PL(A_elecs[i], alpha, E_elecs[i] * Qe, E_max * Qe) * dEe[i] * Qe; //allows comparison of final to initial electron population

        //now calculate the initial critical frequencies
        gamma_e[i] = E_elec_min[i] / Me_EV;  //works as all in eV
        f_em_min[i] = f_crit(gamma_e[i], B); //lower bin values for f_crits

        gamma_e[i] = E_elec_max[i] / Me_EV; //redefine for upper boundary
        f_em_max[i] = f_crit(gamma_e[i], B);

        gamma_e[i] = E_elecs[i] / Me_EV;  //mid points
        beta_e[i] = lortovel(gamma_e[i]); // beta taken at midpts

        f_c[i] = f_crit(gamma_e[i], B);        //assume all radiation occurs at the critical frequency
        dfreqs[i] = f_em_max[i] - f_em_min[i]; //obtain frequency bin widths

        f_pol[i] = f_c[i]; //fixed bin polarization equivalents of the above
        dfreqs_pol[i] = dfreqs[i];

        Ue_jetbase += Ne_e[i] * E_elecs[i] * Qe;
    }

    // SSC Frequency bins
    double bindiff;
    for (i = 0; i < ARRAY_SIZE; i++)
    {
        bindiff = log10(gamma_e[1] * gamma_e[1] * (f_em_min[1])) - log10(gamma_e[0] * gamma_e[0] * (f_em_max[0]));
        f_em_min_IC[i] = pow(10, log10(gamma_e[i] * gamma_e[i] * f_em_min[i]) - bindiff / 2);
        f_em_max_IC[i] = pow(10, log10(gamma_e[i] * gamma_e[i] * f_em_max[i]) + bindiff / 2);
        dfreqs_polIC[i] = (f_em_max_IC[i] - f_em_min_IC[i]);
        f_pol_IC[i] = pow(10, log10(gamma_e[i] * gamma_e[i] * f_pol[i]));
    }

    A_orig = A_elecs[0]; //all elements the same to begin with

    //print the energies to test equipartition
    UB_jetbase = R0 * R0 * C * B0 * B0 / (2 * 1E-7);
    //printf("Obs power, particles, B-fields: %.5e\t%.5e\t%.5e\t%.5e \n", W_j/(gamma_bulk*gamma_bulk), Ue_jetbase, UB_jetbase, Ue_jetbase+UB_jetbase);

    //redefine electron population to ensure equipartition
    for (i = 0; i < ARRAY_SIZE; i++)
    {
        A_elecs[i] = A_PL(alpha, W_j, gamma_bulk, E_min * Qe, E_max * Qe, A_eq) * (UB_jetbase / Ue_jetbase) / A_eq; //1 assumes equipartition. DIVIDED W_j by gamma^2 for jet frame
        dN_dE[i] = electron_PL(A_elecs[i], alpha, E_elecs[i] * Qe, E_max * Qe);
        Ne_e[i] = electron_PL(A_elecs[i], alpha, E_elecs[i] * Qe, E_max * Qe) * dEe[i] * Qe;
        Ne_orig[i] = electron_PL(A_elecs[i], alpha, E_elecs[i] * Qe, E_max * Qe) * dEe[i] * Qe; //allows comparison of final to initial electron population
    }

    //recalaculate equipartition fraction
    Ue_jetbase = 0.0;
    for (i = 0; i < ARRAY_SIZE; i++)
    {
        Ue_jetbase += Ne_e[i] * E_elecs[i] * Qe;
    }
    PRINT("Equipartition fraction U_B / U_e = %.5e \n", UB_jetbase / Ue_jetbase);

    // set bins with fewer than ten electrons to zero: prevents resolution from being limited
    cleanpop(Ne_e, ARRAY_SIZE, 10.0);
    cleanpop(dN_dE, ARRAY_SIZE, 10.0);

    for (i = 0; i < ARRAY_SIZE; i++)
    {
        A_elecs[i] *= (Ne_e[i] / Ne_orig[i]); //sets A to zero for tiny populations
    }

    //11/10/16 define some parameters to determine thick/thin jet sections
    //HELICAL B-field: This variable will save the 3D B-fields on the sky for each jet zone, these B-fields evolve along the jet
    //becoming more transverse
    //defining how far along the helix we are
    size_t helix_counter = 0;
    // sscanf(argv[2], "%d", helix_counter); //taking input from command line how many 'days' we have been observing
    //parameters for helix
    double t_helix = helix_counter * M_PI / 32; //helix parameter x=rcos(t+phi) y=rsin(t+phi) z =ct (r will be radius of jet)
    //Pi constant decided how quickly we sample along the helix -> timescale
    double phi_helix = 0; //M_PI/3; //phase shift, just an initial condition
    double c_helix = R0;  //speed (ie number of coils per unit distance) -- high c => spaced out coils, low c => densely packed coils

    size_t EVPA_rotation = 0;
    // sscanf(argv[1], "%d", EVPA_rotation); //taking input from command line if true begin EVPA rotation, if false do not
    size_t DopDep = 1;
    // sscanf(argv[3], "%d", DopDep); //taking input from command line if true activate DopDep, if false do not

    //now the Bfield vectors are rotated about both the y and x axis each block with a different theta_obs
    //these give the angles each block is at with respect to the direction of the jet:
    //theta_x => rotation angle about x axis, theta_y => rotation angle about y axis
    double th = 2 * atan(tan(deg2rad(theta_open_p)) / gamma_bulk); //full jet diameter opening angle in lab frame in rad
    double theta_r[N_BLOCKS];                                      // Marscher configuration, rings of turbulent cells, middle 3 rings turbulent
    double theta_phi[N_BLOCKS];
    //TEMZ config, calculates theta_r and theta_phi for each block in the concentric circles
    for (i = 0; i < N_BLOCKS; i++)
    {
        if (i < 1)
        {
            theta_r[i] = 0;
            theta_phi[i] = 0;
        }
        else if (i > 0 && i < 7)
        {
            theta_r[i] = th / (2 * (N_RINGS + 0.5));
            theta_phi[i] = -M_PI + (i - 1) * 2 * M_PI / 6;
        }
        else if (i > 6 && i < 19)
        {
            theta_r[i] = 2 * th / (2 * (N_RINGS + 0.5));
            theta_phi[i] = -M_PI + (i - 7) * 2 * M_PI / 12;
        }
        else if (i > 18 && i < 37)
        {
            theta_r[i] = 3 * th / (2 * (N_RINGS + 0.5));
            theta_phi[i] = -M_PI + (i - 19) * 2 * M_PI / 18;
        }
        else if (i > 36 && i < 61)
        {
            theta_r[i] = 4 * th / (2 * (N_RINGS + 0.5));
            theta_phi[i] = -M_PI + (i - 37) * 2 * M_PI / 24;
        }
        else if (i > 60 && i < 91)
        {
            theta_r[i] = 5 * th / (2 * (N_RINGS + 0.5));
            theta_phi[i] = -M_PI + (i - 61) * 2 * M_PI / 30;
        }
        else if (i > 90 && i < 127)
        {
            theta_r[i] = 6 * th / (2 * (N_RINGS + 0.5));
            theta_phi[i] = -M_PI + (i - 91) * 2 * M_PI / 36;
        }
    }

    double theta_tot[N_BLOCKS]; //total angle to line of sight of each block (for doppler boosting) in radians
    for (i = 0; i < N_BLOCKS; i++)
    {
        if (N_BLOCKS == 1)
        {
            theta_tot[i] = deg2rad(theta_obs);
        }
        else
        {
            theta_tot[i] = acos(cos(theta_r[i]) * cos(deg2rad(theta_obs)) + sin(theta_r[i]) * sin(deg2rad(theta_obs)) * cos(M_PI - fabs(theta_phi[i])));
        }
    }

    //alignment vector matrix for each block to every other block, for use in Synchrotron mixing between zones for IC emission, distances are in units of theta_r
    //assuming jet section is flat 2D (ok assumption given large section size) but have z component to rotate with theta_obs
    double align[N_BLOCKS * N_BLOCKS * 3];      //vector from m zone to i zone jet frame
    double unitalign0[N_BLOCKS * N_BLOCKS * 3]; //unit vector from m zone to i zone RPAR adjusted
    double unitalign1[N_BLOCKS * N_BLOCKS * 3]; // unit vector from m zone to i zone RPAR adjusted with offset cosk
    double unitalign2[N_BLOCKS * N_BLOCKS * 3];
    double unitalign3[N_BLOCKS * N_BLOCKS * 3];
    double unitalign4[N_BLOCKS * N_BLOCKS * 3];
    double unitalign5[N_BLOCKS * N_BLOCKS * 3];
    double unitalign6[N_BLOCKS * N_BLOCKS * 3];
    for (i = 0; i < N_BLOCKS; i++)
    {
        for (m = 0; m < N_BLOCKS; m++)
        {
            align[0 + m*3 + i*3*N_BLOCKS] = theta_r[i] * cos(theta_phi[i]) - theta_r[m] * cos(theta_phi[m]); //cos(deg2rad(theta_obs))*(theta_r[i]*cos(theta_phi[i]) - theta_r[m]*cos(theta_phi[m]));
            align[1 + m*3 + i*3*N_BLOCKS] = theta_r[i] * sin(theta_phi[i]) - theta_r[m] * sin(theta_phi[m]);
            align[2 + m*3 + i*3*N_BLOCKS] = 0.0; //-sin(deg2rad(theta_obs))*(theta_r[i]*cos(theta_phi[i]) - theta_r[m]*cos(theta_phi[m]));
        }
    }

    size_t breakstep = 100;
    // sscanf(argv[8], "%d", breakstep);
    //Full and zonal Stokes vector bins for Synchrotron and IC respectively
    double S_Stokes[N_BLOCKS][3][ARRAY_SIZE];
    memset(S_Stokes, 0, sizeof(S_Stokes[0][0][0]) * N_BLOCKS * ARRAY_SIZE * 3);
    double S_StokesTotal_Sec[breakstep][ARRAY_SIZE][3];
    memset(S_StokesTotal_Sec, 0, sizeof(S_StokesTotal_Sec[0][0][0]) * breakstep * ARRAY_SIZE * 3);

    double zone_doppler;
    double u_rad = 0.0;
    double Blength; //for normalization
    double logdif;
    int binshiftIC;
    int binshiftS;
    double dfactor;

    double q_theta_sync;
    double opacity_factor;
    double dx_op_array[breakstep];
    double k_array[breakstep][ARRAY_SIZE];
    double tau[breakstep][ARRAY_SIZE];
    memset(tau, 0, sizeof(tau[0][0]) * breakstep * ARRAY_SIZE);

    //choosing random vectors (B_0,B_1,B_2) in unit sphere for random blocks initial B directions:
    MTRand seedr = seedRand(argblocks[2]); //seedRand((unsigned)time(NULL)+(unsigned)(1631*task_id + 335*thread_id)); //random seed supplemented by task and thread id
    //MTRand seedr = seedRand(11); //fix random seed
    
    double Bx_helical[N_BLOCKS]; //simplifies angular shifts for the helical components
    double By_helical[N_BLOCKS];
    double Bz_helical[N_BLOCKS];
    double BX[nLT_sections][N_BLOCKS]; //for joint helical and random array
    double BY[nLT_sections][N_BLOCKS];
    double BZ[nLT_sections][N_BLOCKS];
    //B_effectives due to RPAR for each zone in every other zone
    double B_effectives[nLT_sections * N_BLOCKS * N_BLOCKS * 3];
    memset(B_effectives, 0, sizeof(B_effectives[0]) * nLT_sections * N_BLOCKS * N_BLOCKS * 3);

    int marker_list[N_BLOCKS];
    for (n = 0; n < N_BLOCKS; n++)
    {
        marker_list[n] = nLT_sections / 2;
    }

    double B_0[nLT_sections][N_BLOCKS];
    double B_1[nLT_sections][N_BLOCKS];
    double B_2[nLT_sections][N_BLOCKS];
    double unif_theta;
    double unif_phi;
    for (i = 0; i < (N_BLOCKS); i++)
    { //these are the random B-field vectors in each block as if we are looking straight down jet (have to rotate by theta_obs)
        for (l = 0; l < nLT_sections; l++)
        {
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
    for (i = 0; i < (N_BLOCKS); i++)
    {
        Bx_helical[i] = -sin(t_helix + phi_helix) / sqrt(2);
        By_helical[i] = cos(t_helix + phi_helix) / sqrt(2); //no transversification or cone surface, 45 deg helix
        Bz_helical[i] = 1 / sqrt(2);
    }

    if (!EVPA_rotation)
    {
        for (l = 0; l < nLT_sections; l++)
        {
            for (i = 0; i < (N_BLOCKS); i++)
            {
                BX[l][i] = (B_2[l][i] * sin(deg2rad(theta_obs)) + B_0[l][i] * cos(deg2rad(theta_obs)));
                BY[l][i] = B_1[l][i];
                BZ[l][i] = (B_2[l][i] * cos(deg2rad(theta_obs)) - B_0[l][i] * sin(deg2rad(theta_obs)));
            }
        }
    }
    else
    {
        for (i = 0; i < (N_BLOCKS); i++)
        {
            if (i < 19)
            {
                BX[l][i] = (Bz_helical[i] * sin(deg2rad(theta_obs)) + Bx_helical[i] * cos(deg2rad(theta_obs)));
                BY[l][i] = By_helical[i];
                BZ[l][i] = (Bz_helical[i] * cos(deg2rad(theta_obs)) - Bx_helical[i] * sin(deg2rad(theta_obs)));
            }
            else
            {
                BX[l][i] = (B_2[l][i] * sin(deg2rad(theta_obs)) + B_0[l][i] * cos(deg2rad(theta_obs)));
                BY[l][i] = B_1[l][i];
                BZ[l][i] = (B_2[l][i] * cos(deg2rad(theta_obs)) - B_0[l][i] * sin(deg2rad(theta_obs)));
            }
        }
    }

    if (!DopDep)
    {
        for (i = 0; i < (N_BLOCKS); i++)
        {
            for (l = 0; l < N_BLOCKS; l++)
            {
                Blength = sqrt(BX[m][i] * BX[m][i] + BY[m][i] * BY[m][i] + BZ[m][i] * BZ[m][i]); //normalization important to get correct z angle

                B_effectives[0 + i*3 + l*3*N_BLOCKS + m*3*N_BLOCKS*N_BLOCKS] = BX[m][i] / Blength;
                B_effectives[1 + i*3 + l*3*N_BLOCKS + m*3*N_BLOCKS*N_BLOCKS] = BY[m][i] / Blength;
                B_effectives[2 + i*3 + l*3*N_BLOCKS + m*3*N_BLOCKS*N_BLOCKS] = BZ[m][i] / Blength;
            }
        }
    }
    else
    {
        for (i = 0; i < N_BLOCKS; i++)
        {
            for (l = 0; l < (N_BLOCKS); l++)
            {

                for (m = 0; m < nLT_sections; m++)
                {
                    DD_3Beffective(BX[m][i], BY[m][i], BZ[m][i],
                                   cos(deg2rad(theta_obs)) * sin(theta_r[l]) * cos(theta_phi[l]) + sin(deg2rad(theta_obs)) * cos(theta_r[l]),
                                   sin(theta_r[l]) * sin(theta_phi[l]), -sin(deg2rad(theta_obs)) * sin(theta_r[l]) * cos(theta_phi[l]) + cos(deg2rad(theta_obs)) * cos(theta_r[l]),
                                   gamma_bulk, &B_effectives[i*3 + l*3*N_BLOCKS + m*3*N_BLOCKS*N_BLOCKS]);
                }

                DD_3Beffective(align[0 + i*3 + l*3*N_BLOCKS], align[1 + i*3 + l*3*N_BLOCKS], align[2 + i*3 + l*3*N_BLOCKS],
                               cos(deg2rad(theta_obs)) * sin(theta_r[l]) * cos(theta_phi[l]) + sin(deg2rad(theta_obs)) * cos(theta_r[l]),
                               sin(theta_r[l]) * sin(theta_phi[l]), -sin(deg2rad(theta_obs)) * sin(theta_r[l]) * cos(theta_phi[l]) + cos(deg2rad(theta_obs)) * cos(theta_r[l]),
                               gamma_bulk, &unitalign0[i*3 + l*3*N_BLOCKS]);

                DD_3Beffective(align[0 + i*3 + l*3*N_BLOCKS], align[1 + i*3 + l*3*N_BLOCKS], sin(18.2 * M_PI / 180) * sqrt(pow(align[0 + i*3 + l*3*N_BLOCKS], 2) + pow(align[1 + i*3 + l*3*N_BLOCKS], 2)), //adding cosk offsets
                               cos(deg2rad(theta_obs)) * sin(theta_r[l]) * cos(theta_phi[l]) + sin(deg2rad(theta_obs)) * cos(theta_r[l]),
                               sin(theta_r[l]) * sin(theta_phi[l]), -sin(deg2rad(theta_obs)) * sin(theta_r[l]) * cos(theta_phi[l]) + cos(deg2rad(theta_obs)) * cos(theta_r[l]),
                               gamma_bulk, &unitalign1[i*3 + l*3*N_BLOCKS]);

                DD_3Beffective(align[0 + i*3 + l*3*N_BLOCKS], align[1 + i*3 + l*3*N_BLOCKS], sin(36.4 * M_PI / 180) * sqrt(pow(align[0 + i*3 + l*3*N_BLOCKS], 2) + pow(align[1 + i*3 + l*3*N_BLOCKS], 2)),
                               cos(deg2rad(theta_obs)) * sin(theta_r[l]) * cos(theta_phi[l]) + sin(deg2rad(theta_obs)) * cos(theta_r[l]),
                               sin(theta_r[l]) * sin(theta_phi[l]), -sin(deg2rad(theta_obs)) * sin(theta_r[l]) * cos(theta_phi[l]) + cos(deg2rad(theta_obs)) * cos(theta_r[l]),
                               gamma_bulk, &unitalign2[i*3 + l*3*N_BLOCKS]);

                DD_3Beffective(align[0 + i*3 + l*3*N_BLOCKS], align[1 + i*3 + l*3*N_BLOCKS], sin(-18.2 * M_PI / 180) * sqrt(pow(align[0 + i*3 + l*3*N_BLOCKS], 2) + pow(align[1 + i*3 + l*3*N_BLOCKS], 2)),
                               cos(deg2rad(theta_obs)) * sin(theta_r[l]) * cos(theta_phi[l]) + sin(deg2rad(theta_obs)) * cos(theta_r[l]),
                               sin(theta_r[l]) * sin(theta_phi[l]), -sin(deg2rad(theta_obs)) * sin(theta_r[l]) * cos(theta_phi[l]) + cos(deg2rad(theta_obs)) * cos(theta_r[l]),
                               gamma_bulk, &unitalign3[i*3 + l*3*N_BLOCKS]);

                DD_3Beffective(align[0 + i*3 + l*3*N_BLOCKS], align[1 + i*3 + l*3*N_BLOCKS], sin(-36.4 * M_PI / 180) * sqrt(pow(align[0 + i*3 + l*3*N_BLOCKS], 2) + pow(align[1 + i*3 + l*3*N_BLOCKS], 2)),
                               cos(deg2rad(theta_obs)) * sin(theta_r[l]) * cos(theta_phi[l]) + sin(deg2rad(theta_obs)) * cos(theta_r[l]),
                               sin(theta_r[l]) * sin(theta_phi[l]), -sin(deg2rad(theta_obs)) * sin(theta_r[l]) * cos(theta_phi[l]) + cos(deg2rad(theta_obs)) * cos(theta_r[l]),
                               gamma_bulk, &unitalign4[i*3 + l*3*N_BLOCKS]);

                DD_3Beffective(align[0 + i*3 + l*3*N_BLOCKS], align[1 + i*3 + l*3*N_BLOCKS], sin(-60 * M_PI / 180) * sqrt(pow(align[0 + i*3 + l*3*N_BLOCKS], 2) + pow(align[1 + i*3 + l*3*N_BLOCKS], 2)),
                               cos(deg2rad(theta_obs)) * sin(theta_r[l]) * cos(theta_phi[l]) + sin(deg2rad(theta_obs)) * cos(theta_r[l]),
                               sin(theta_r[l]) * sin(theta_phi[l]), -sin(deg2rad(theta_obs)) * sin(theta_r[l]) * cos(theta_phi[l]) + cos(deg2rad(theta_obs)) * cos(theta_r[l]),
                               gamma_bulk, &unitalign5[i*3 + l*3*N_BLOCKS]);

                DD_3Beffective(align[0 + i*3 + l*3*N_BLOCKS], align[1 + i*3 + l*3*N_BLOCKS], sin(60 * M_PI / 180) * sqrt(pow(align[0 + i*3 + l*3*N_BLOCKS], 2) + pow(align[1 + i*3 + l*3*N_BLOCKS], 2)),
                               cos(deg2rad(theta_obs)) * sin(theta_r[l]) * cos(theta_phi[l]) + sin(deg2rad(theta_obs)) * cos(theta_r[l]),
                               sin(theta_r[l]) * sin(theta_phi[l]), -sin(deg2rad(theta_obs)) * sin(theta_r[l]) * cos(theta_phi[l]) + cos(deg2rad(theta_obs)) * cos(theta_r[l]),
                               gamma_bulk, &unitalign6[i*3 + l*3*N_BLOCKS]);
            }
        }
    }

    //Linear Algebra Speedup using BLAS and GSL libraries: write SSC integral as y = Ax, code below finds A

    double effective_alpha[ARRAY_SIZE];
    calculate_alpha(dN_dE, E_elecs, effective_alpha);
    
    gsl_vector *X[N_BLOCKS];
    gsl_matrix *A[N_BLOCKS*3];
    gsl_vector *IC_Stokes[N_BLOCKS*3];

    for (h = 0; h < N_BLOCKS; h++)
    {
        X[h] = gsl_vector_calloc(ARRAY_SIZE * ARRAY_SIZE * N_BLOCKS * 2); //perp and para concatenated

        A[0+h*3] = gsl_matrix_calloc(ARRAY_SIZE, ARRAY_SIZE * ARRAY_SIZE * N_BLOCKS * 2); //one for each stokes parameter
        A[1+h*3] = gsl_matrix_calloc(ARRAY_SIZE, ARRAY_SIZE * ARRAY_SIZE * N_BLOCKS * 2);
        A[2+h*3] = gsl_matrix_calloc(ARRAY_SIZE, ARRAY_SIZE * ARRAY_SIZE * N_BLOCKS * 2);

        IC_Stokes[0+3*h] = gsl_vector_calloc(ARRAY_SIZE);
        IC_Stokes[1+3*h] = gsl_vector_calloc(ARRAY_SIZE);
        IC_Stokes[2+3*h] = gsl_vector_calloc(ARRAY_SIZE);
    }
    
    // PRINT("Marker list: \n", nLT_sections/2);
    if (SSC)
    {   
        set_SSC_matrix(A, N_BLOCKS, f_pol, f_pol_IC, 
                    B_effectives,
                    E_elecs, E_elec_min, E_elec_max,
                    dEe, effective_alpha,
                    unitalign0, unitalign1, unitalign2,
                    unitalign3, unitalign4, unitalign5,
                    unitalign6);
    }

    // Actual Jet emission calculation loop
    // printf("Beginning jet analysis... \n");
    while (x < L_jet)
    {
        eps = epsilon(B); //same for all populations

        //----------------- B-field projection onto sky for this section --------------//
        // Now have introduced R dependence on the random blocks as well, works just like for the helical case, just * R_old/R_new on B_2 term
        //now can also include doppler depolarisation, changing the B field to EFFECTIVE B fields which give correct perp E component as if it had
        //been doppler depolarisedx
        //Each block has its own theta_obs as we are looking inside the conical jet

        UB_jetbase = R * R * C * B * B / (2 * 1E-7);
        Ue_jetbase = 0.0;
        for (i = 0; i < ARRAY_SIZE; i++)
        {
            Ue_jetbase += Ne_e[i] * E_elecs[i] * Qe;
        }
        // printf("Ue = %.5e\n",Ue_jetbase);
        // printf("UB = %.5e\n",UB_jetbase);
        // printf("Equipartition fraction U_B / U_e = %.3e\n", UB_jetbase/Ue_jetbase);

        //Jet section mixing here: if Rtan(theta_obs) > R_0/sqrt(N) then we start assigning new B-fields to edge zones and so on, equivalent to seeing emission from zones in section
        //behind and in front.
        if (x > 2 * R0 / (sqrt(N_BLOCKS)) / (1 - beta_bulk))
        { //test to make sure mixing due to finite gamma doesnt occur
            printf("\n WARNING: GAMMA MIXING! \n");
        }

        //******************************************************************************************************//
        //Light travel angular effect + gamma bulk effect: causes section mixing
        if (R * tan(deg2rad(theta_obs)) > R0 / (sqrt(N_BLOCKS)))
        {
            for (i = 0; i < (N_BLOCKS); i++)
            {
                if (fabs(R * 2 * theta_r[i] * cos(theta_phi[i]) / th) > (abs(2 * (marker_list[i] - 14) - 1) * R0 / (sqrt(N_BLOCKS))) / tan(deg2rad(theta_obs)))
                {
                    if (fabs(theta_phi[i]) > M_PI / 2)
                    {
                        marker_list[i] += 1;
                    }
                    else
                    {
                        if (marker_list[i] == 15)
                        {
                            marker_list[i] -= 2;
                        }
                        else
                        {
                            marker_list[i] -= 1;
                        }
                    }
                }
                //printf("i, zone_marker = %d \t  %d \n", i, marker_list[i]);
            }
            // printf("\n WARNING: Section Mixing has begun! \n");
        }
        /***************************************************************************************************/

        //begin with Synchrotron emission
        for (i = 0; i < ARRAY_SIZE; i++)
        {
            P_single[i] = larmor(B, beta_e[i], gamma_e[i]); //gives the radiated power by one electron in each bin
            //j[i] = j_per_hz(P_single[i], R, dN_dE[i], eps, f_c[i]);//*(R_next/R)*(R_next/R);//gives emissivity per unit frequency

            //k[i] = k_new(j[i], eps, f_c[i]); //gives the corresponding opacity
            Ps_per_m_elec[i] = P_single[i] * dN_dE[i] * E_elecs[i] * Qe * dfreqs_pol[i] / (2 * f_c[i] * C);
            //Ps_per_m_test[i] = power_emitted(j[i], k[i], R) * M_PI * R; //power per unit length assuming R roughly constant
            //determine whether jet is optically thin or thick -> this part is incomplete
            //            R_eff[i] = power_emitted(j[i], k[i], R)/j[i]; //depth down to which can be seen in one second

            //Polarisation -> Power emitted parallel and perpendicular to projected B field at each frequency interval per m (/ l_c)
            for (l = 0; l < ARRAY_SIZE; l++)
            {
                for (m = 0; m < 75; m++)
                {
                    if ((f_pol[l] / f_c[i]) >= FG[m][0] && (f_pol[l] / f_c[i]) <= FG[m + 1][0])
                    {
                        P_perp[l] += (sqrt(3) * M_PI / (16 * Me * pow(C, 3))) * pow(Qe, 3) * B * FG[m + 1][1] * dN_dE[i] * dEe[i] * Qe * /*dfreqs_pol[l] */ 4.4216E11 / (1E-7); //power per unit frequency, this is Power(freq), multiply by freq later to get nuF(nu)
                        P_para[l] += (sqrt(3) * M_PI / (16 * Me * pow(C, 3))) * pow(Qe, 3) * B * FG[m + 1][2] * dN_dE[i] * dEe[i] * Qe * /*dfreqs_pol[l] */ 4.4216E11 / (1E-7); //*dx !!! (this happens further down) and * other constants, already have /l_c
                    }

                    else if ((f_pol[l] / f_c[i]) <= FG[0][0] || (f_pol[l] / f_c[i]) >= FG[73][0])
                    {
                        P_perp[l] += 0.0;
                        P_para[l] += 0.0;
                    }
                }
            }
            Sync_losses[i] = Ps_per_m_elec[i];
        }

        for (i = 0; i < ARRAY_SIZE; i++)
        {
            j[i] = (P_perp[i] + P_para[i]) / (M_PI * R * R);
            k[i] = (j[i] * C * C) / (2 * pow(eps, 0.5) * pow(f_pol[i], 2.5));
        }


        //some code to estimate what dx early for energy density calcultion, based off sync losses only, ok if IC losses not too big, early jet
        dx_set = findminelementNO0(Ne_e, ARRAY_SIZE); //gets the lowest non_zero element. Higher energy electrons radiate more rapidly

        dx_R = 0.05 * (R0 + x * tan(deg2rad(theta_open_p))) / tan(deg2rad(theta_open_p));               //ensures Rnew <= 1.05 Rold
        dx_P = (Ne_e[dx_set] * (E_elecs[dx_set] - E_elecs[dx_set - 1]) * Qe) / (Ps_per_m_elec[dx_set]); //based on radiative losses

        if (dx_R < dx_P)
        {
            dx = dx_R;
        }
        else
        {
            dx = dx_P;
        }
        //
        //        for (i=0; i<ARRAY_SIZE; i++){
        //            printf("Opacities %.5e\t%.5e\n", k[i], (1 - exp(-k[i]*dx)) / k[i]);
        //        }

        /*************************************************************************************************************/
        // Buffer of previous synchrotron emission setup
        if (x == 0)
        {
            // max buffer_size is when x == 2R_array[min]
            R_array = (double *)malloc(sizeof(double));
            dx_array = (double *)malloc(sizeof(double));
            Pperp_array = (double **)malloc(1 * sizeof(double *));
            Pperp_array[0] = (double *)malloc(ARRAY_SIZE * sizeof(double)); //dynamically allocating 2D double arrays
            Ppara_array = (double **)malloc(1 * sizeof(double *));
            Ppara_array[0] = (double *)malloc(ARRAY_SIZE * sizeof(double));

            if (R_array == NULL || dx_array == NULL)
            {
                printf("malloc failed\n");
                fprintf(stderr, "malloc failed\n");
                return (-1);
            }

            R_array[0] = R0;
            dx_array[0] = dx;
            for (l = 0; l < ARRAY_SIZE; l++)
            {
                Ppara_array[0][l] = P_para[l];
                Pperp_array[0][l] = P_perp[l];
            }
        }
        else if (2 * theta_r[N_BLOCKS - 1] / th * R_array[0] + R_array[findminelement(R_array, buffer_size)] > x + dx)
        { //add rows to matrices, else do nothing
            for (n = 0; n < N_BLOCKS; n++)
            {
                if (2 * theta_r[n] / th * R_array[0] + R_array[findminelement(R_array, buffer_size)] > x + dx)
                {
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
        if (x != 0 && buffer_size > 1)
        {
            for (n = buffer_size - 2; n >= 0; n--)
            {
                R_array[n + 1] = R_array[n];
                dx_array[n + 1] = dx_array[n];
                for (l = 0; l < ARRAY_SIZE; l++)
                {   
                    Pperp_array[n + 1][l] = Pperp_array[n][l];
                    Ppara_array[n + 1][l] = Ppara_array[n][l];
                }
            }

            R_array[0] = R;
            dx_array[0] = dx;
            for (l = 0; l < ARRAY_SIZE; l++)
            {
                Pperp_array[0][l] = P_perp[l];
                Ppara_array[0][l] = P_para[l];
            }
        }

        /******************************************************************************************************************/
        Area0 = circle_area(2 * theta_r[0] / th, dx_array[0], dx_array[0], R_array[0]);

        // Calculate effective alpha for electron population as function of energy and Calculate Synchrotron Stokes parameters

        calculate_alpha(dN_dE, E_elecs, effective_alpha);

        for (l = 0; l < ARRAY_SIZE; l++)
        {
            //printf("alphaef \t%.5e\n",effective_alpha[l]);
            for (g = 0; g < N_BLOCKS; g++)
            {   //Sync Stokes Parameters
                q_theta_sync = pow(sin(acos(B_effectives[2 + g*3 + g*3*N_BLOCKS + marker_list[g]*3*N_BLOCKS*N_BLOCKS])), (effective_alpha[l] + 1) / 2); //factor to account for weaker emission when B_field pointed closer to line of sight
                S_Stokes[g][0][l] += q_theta_sync * P_perp[l] / N_BLOCKS;
                S_Stokes[g][1][l] += (q_theta_sync * P_perp[l] / N_BLOCKS) * cos(2 * atan2(B_effectives[0 + g*3 + g*3*N_BLOCKS + marker_list[g]*3*N_BLOCKS*N_BLOCKS], -B_effectives[1 + g*3 + g*3*N_BLOCKS + marker_list[g]*3*N_BLOCKS*N_BLOCKS]));
                S_Stokes[g][2][l] += (q_theta_sync * P_perp[l] / N_BLOCKS) * sin(2 * atan2(B_effectives[0 + g*3 + g*3*N_BLOCKS + marker_list[g]*3*N_BLOCKS*N_BLOCKS], -B_effectives[1 + g*3 + g*3*N_BLOCKS + marker_list[g]*3*N_BLOCKS*N_BLOCKS]));
                S_Stokes[g][0][l] += q_theta_sync * P_para[l] / N_BLOCKS;
                S_Stokes[g][1][l] += (q_theta_sync * P_para[l] / N_BLOCKS) * cos(2 * atan2(B_effectives[1 + g*3 + g*3*N_BLOCKS + marker_list[g]*3*N_BLOCKS*N_BLOCKS], B_effectives[0 + g*3 + g*3*N_BLOCKS + marker_list[g]*3*N_BLOCKS*N_BLOCKS]));
                S_Stokes[g][2][l] += (q_theta_sync * P_para[l] / N_BLOCKS) * sin(2 * atan2(B_effectives[1 + g*3 + g*3*N_BLOCKS + marker_list[g]*3*N_BLOCKS*N_BLOCKS], B_effectives[0 + g*3 + g*3*N_BLOCKS + marker_list[g]*3*N_BLOCKS*N_BLOCKS]));
            }

            //Energy Density loop, calculates the total energy density seen at each zone for this section of the jet
            // Always need to have uradARRAY[] = sum (dfactor over second index)
            // dfactor array should be all that is needed, dfactor[g][h][l] takes over role of urad_array_perp[l]
            // dfactor[g][h][l] = synchrotron energy density at zone h from zone g at frequency l
            for (n = 0; n < N_BLOCKS; n++)
            {
                dxsum = 0;
                for (i = 0; i < buffer_subset[n]; i++)
                {
                    if (dxsum + dx_array[i] < 2 * theta_r[n] / th * R_array[0] + R_array[i])
                    { //this takes place of buffer_subset, making sure emission from further than R_i not contributing
                        dxsum += dx_array[i];
                        urad_array_perp[n][l] += (Pperp_array[i][l]) * dfreqs_pol[l] * dx_array[i] / (M_PI * pow(R_array[i], 2) * C) * (circle_area(2 * theta_r[n] / th * R_array[0], dx_array[i], dxsum, R_array[i]) / Area0) * (dx_array[0] / dxsum);
                        urad_array_para[n][l] += (Ppara_array[i][l]) * dfreqs_pol[l] * dx_array[i] / (M_PI * pow(R_array[i], 2) * C) * (circle_area(2 * theta_r[n] / th * R_array[0], dx_array[i], dxsum, R_array[i]) / Area0) * (dx_array[0] / dxsum);

                        counter = 0.0; //divides up number zones contributing at the same time
                        for (g = 0; g < N_BLOCKS; g++)
                        {
                            dist = fabs(sqrt(pow(2 * align[g*3*N_BLOCKS] / th * R_array[i] - 2 * theta_r[n] / th * R_array[0] * cos(theta_phi[n]), 2) + pow(2 * align[1 + g*3*N_BLOCKS] / th * R_array[i] - 2 * theta_r[n] / th * R_array[0] * sin(theta_phi[n]), 2)) - dxsum);
                            if (dist < 1 / (2 * (N_RINGS + 0.5)) * R_array[i])
                            {
                                counter += 1;
                            }
                        }

                        for (g = 0; g < N_BLOCKS; g++)
                        {
                            dist = fabs(sqrt(pow(2 * align[g*3*N_BLOCKS] / th * R_array[i] - 2 * theta_r[n] / th * R_array[0] * cos(theta_phi[n]), 2) + pow(2 * align[1 + g*3*N_BLOCKS] / th * R_array[i] - 2 * theta_r[n] / th * R_array[0] * sin(theta_phi[n]), 2)) - dxsum);
                            if (dist < 1 / (2 * (N_RINGS + 0.5)) * R_array[i])
                            {
                                dfactor_perp[n][g][l] += (Pperp_array[i][l]) * dfreqs_pol[l] * dx_array[i] / (M_PI * pow(R_array[i], 2) * C) * (circle_area(2 * theta_r[n] / th * R_array[0], dx_array[i], dxsum, R_array[i]) / Area0) * (dx_array[0] / dxsum) / counter;
                                dfactor_para[n][g][l] += (Ppara_array[i][l]) * dfreqs_pol[l] * dx_array[i] / (M_PI * pow(R_array[i], 2) * C) * (circle_area(2 * theta_r[n] / th * R_array[0], dx_array[i], dxsum, R_array[i]) / Area0) * (dx_array[0] / dxsum) / counter;
                            }
                        }
                    }
                }
            }
        }

        //********************************************** ALP SSC losses *******************************************//
        // given energy density prescription do we need 1/nblocks in Pperpperp etc? Yes, still do because Pperp_array uses full electron pop power
        if (SSC)
        {
            for (h = 0; h < N_BLOCKS; h++)
            { //B-fields
                for (g = 0; g < N_BLOCKS; g++)
                { //B-fields
                    for (l = 0; l < ARRAY_SIZE; l++)
                    { //Seed Frequencies
                        for (i = 0; i < ARRAY_SIZE; i++)
                        { //Scattering electron energy bins
                            gsl_vector_set(X[h], (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2 + 1, dN_dE[i] * dfactor_perp[h][g][l]);
                            gsl_vector_set(X[h], (i + ARRAY_SIZE * l + ARRAY_SIZE * ARRAY_SIZE * g) * 2, dN_dE[i] * dfactor_para[h][g][l]);
                        }
                    }
                }
            }
            // printf("dfactor: %f\t%f\t%f\n", log10(dfactor_perp[0][0][10]), log10(dfactor_perp[0][0][25]), log10(dfactor_perp[0][0][40]));
            // printf("dNdE: %f\t%f\t%f\n", log10(dN_dE[10]), log10(dN_dE[25]), log10(dN_dE[40]));

            for (h = 0; h < N_BLOCKS; h++)
            {
                gsl_blas_dgemv(CblasNoTrans, 1, A[0+h*3], X[h], 0, IC_Stokes[0+h*3]);
                gsl_blas_dgemv(CblasNoTrans, 1, A[1+h*3], X[h], 0, IC_Stokes[1+h*3]);
                gsl_blas_dgemv(CblasNoTrans, 1, A[2+h*3], X[h], 0, IC_Stokes[2+h*3]);
            }
        }

        for (i = 0; i < ARRAY_SIZE; i++)
        {
            for (n = 0; n < N_BLOCKS; n++)
            {
                u_rad += (urad_array_perp[n][i] + urad_array_para[n][i]);
            }
        }
        u_rad /= N_BLOCKS; //use average energy density for electron losses, otherwise different blocks cool at different rates - SSC cooling isnt that important anyway
        // PRINT("u_rad %.5e \n", u_rad);
        for (i = 0; i < ARRAY_SIZE; i++) //1.79 //0.85
        {
            Ps_per_mIC_elec[i] = 0.0;
            Ps_per_mIC_elec[i] = (4 / 3) * 6.6524E-29 * gamma_e[i] * gamma_e[i] * beta_e[i] * beta_e[i] * dN_dE[i] * dEe[i] * Qe * u_rad; //test for Ps_per_mIC_elec == this
            IC_losses[i] = Ps_per_mIC_elec[i];
        }

        //**********************************************************************************************//
        //some code to estimate what dx should be
        dx_set = findminelementNO0(Ne_e, ARRAY_SIZE); //gets the lowest non_zero element. Higher energy electrons radiate more rapidly

        dx_R = 0.05 * (R0 + x * tan(deg2rad(theta_open_p))) / tan(deg2rad(theta_open_p));                                         //ensures Rnew <= 1.05 Rold
        dx_P = (Ne_e[dx_set] * (E_elecs[dx_set] - E_elecs[dx_set - 1]) * Qe) / (Ps_per_m_elec[dx_set] + Ps_per_mIC_elec[dx_set]); //based on radiative losses

        if (dx_R < dx_P)
        {
            dx = dx_R;
        }
        else
        {
            dx = dx_P;
        }

        for (h = 0; h < N_BLOCKS; h++)
        { //summing the Stokes parameters of all the blocks and boosting
            zone_doppler = doppler(beta_bulk, theta_tot[h]);
            //difference between bin boundaries in logspace:
            logdif = log10(f_pol_IC[1]) - log10(f_pol_IC[0]);
            //adjust zone_doppler with magnetic field z value as in Marscher
            double adjust = 1.;
            if (argblocks[3] == 1)
            {
                adjust = pow(B_Z[marker_list[h]][h],2) / 0.33333;
            }
            binshiftIC = (int)round(log10(zone_doppler / doppler_factor) / logdif); //the number of bins one moves across (relative to shift of bins which already takes place)
            logdif = log10(f_pol[1]) - log10(f_pol[0]);
            binshiftS = (int)round(log10(zone_doppler * adjust / doppler_factor) / logdif);
            // printf("binshifts %d\t%d\n", binshiftIC,binshiftS);
            // printf("dopplers %.5e\t%.5e\t%.5e\n", zone_doppler,doppler_factor,logdif);
            // printf("IC_Stokes %.5e\t%.5e\t%.5e\n", gsl_vector_get(IC_Stokes[0],19),
            //                                         gsl_vector_get(IC_Stokes[1],19),
            //                                         gsl_vector_get(IC_Stokes[2],19));

            for (n = 0; n < ARRAY_SIZE; n++)
            {
                // Might want zone_doppler ^ 3 for a continuous jet
                if (n - binshiftIC >= 0 && n - binshiftIC < ARRAY_SIZE)
                {
                    IC_StokesTotal[n + 0 * ARRAY_SIZE] += gsl_vector_get(IC_Stokes[0+h*3], n - binshiftIC) * pow(zone_doppler, 4) * dx;
                    IC_StokesTotal[n + 1 * ARRAY_SIZE] += gsl_vector_get(IC_Stokes[1+h*3], n - binshiftIC) * pow(zone_doppler, 4) * dx;
                    IC_StokesTotal[n + 2 * ARRAY_SIZE] += gsl_vector_get(IC_Stokes[2+h*3], n - binshiftIC) * pow(zone_doppler, 4) * dx;
                }
                if (n - binshiftS >= 0 && n - binshiftS < ARRAY_SIZE)
                {
                    opacity_factor = (1 - exp(-k[n - binshiftS] * dx)) / k[n - binshiftS]; //tends to dx as opacity goes to zero (optically thin)
                    if (opacity_factor > 0.9 * dx || isnan(opacity_factor) || opacity_factor == 0.0)
                    {
                        opacity_factor = dx;
                    }

                    S_StokesTotal_Sec[nSteps][n][0] += S_Stokes[h][0][n - binshiftS] * pow(zone_doppler, 4) * opacity_factor; //* dx;
                    S_StokesTotal_Sec[nSteps][n][1] += S_Stokes[h][1][n - binshiftS] * pow(zone_doppler, 4) * opacity_factor; //* dx;
                    S_StokesTotal_Sec[nSteps][n][2] += S_Stokes[h][2][n - binshiftS] * pow(zone_doppler, 4) * opacity_factor; //* dx;
                }
            }
        }
        //change values for next loop
        x += dx;
        R_prev = R;
        R = R_new(R0, deg2rad(theta_open_p), x);

        //correct synchrotron for full length + IC
        for (i = 0; i < ARRAY_SIZE; i++)
        {
            Ps_per_m[i] *= dx;
            //Ps_per_m_IC[i] *= dx;
            opacity_factor = (1 - exp(-k[i] * dx)) / k[i]; //tends to dx as opacity goes to zero (optically thin)
            if (opacity_factor > 0.9 * dx || isnan(opacity_factor) || opacity_factor == 0.0)
            {
                opacity_factor = dx;
            }
            Sync_losses[i] *= opacity_factor; //low energy synchrotron losses now affected by opacity
            IC_losses[i] *= dx;               //ALP
            P_perp[i] *= dx;                  //polarisation powers
            P_para[i] *= dx;
        }

        //compute the synchrotron losses + IC LOSSES
        for (i = 0; i < ARRAY_SIZE; i++)
        {
            if (i == 0)
            {
                //Ne_losses[i]=0.0; //cannot drop to lower bins
                Ne_losses[i] = (Sync_losses[i] + IC_losses[i]) / ((E_elecs[i] - Me_EV) * Qe); //can't radiate anymore
            }
            else
            {
                // Ne_losses[i] = (Sync_losses[i])*(C/C)/((E_elecs[i]-E_elecs[i-1])*Qe);
                Ne_losses[i] = (Sync_losses[i] + IC_losses[i]) / ((E_elecs[i] - E_elecs[i - 1]) * Qe);
            }
            if (i == ARRAY_SIZE - 1)
            {
                Ne_gains[i] = 0.0;
            }
            else
            {
                //Ne_gains[i] = (Sync_losses[i+1])*(C/C)/((E_elecs[i+1]-E_elecs[i])*Qe);
                Ne_gains[i] = (Sync_losses[i + 1] + IC_losses[i + 1]) * (C / C) / ((E_elecs[i + 1] - E_elecs[i]) * Qe);
            }
        }

        for (i = 0; i < ARRAY_SIZE; i++)
        {
            Ne_e[i] = Ne_e[i] - Ne_losses[i] + Ne_gains[i];
        }

        cleanpop(Ne_e, ARRAY_SIZE, 10.0);        //any elements with Ne_e<10 set to zero to avoid nans

        B_prev = B; //need to update dfreqs
        B = get_newB(B0, R0, R);
        eps = epsilon(B);

        // reset values to 0 for the next section x step
        u_rad = 0.0;
        for (i = 0; i < ARRAY_SIZE; i++)
        {
            f_c[i] = f_crit(gamma_e[i], B);
            f_em_min[i] *= B / B_prev; //steps[nSteps];
            f_em_max[i] *= B / B_prev; //steps[nSteps];
            dfreqs[i] = f_em_max[i] - f_em_min[i];
            A_elecs[i] = A_orig * (Ne_e[i] / Ne_orig[i]);
            dN_dE[i] = electron_PL(A_elecs[i], alpha, E_elecs[i] * Qe, E_max * Qe);
            k_array[nSteps][i] = k[i]; //opacity of section for later integral
            P_para[i] = 0.0;           //Reset
            P_perp[i] = 0.0;
            Ps_per_mIC_elec[i] = 0.0;
            Ps_per_m_elec[i] = 0.0;
        }

        memset(S_Stokes, 0, sizeof(S_Stokes[0][0][0]) * N_BLOCKS * ARRAY_SIZE * 3); //reset these to 0 for next section
        //memset(IC_Stokes, 0, sizeof(IC_Stokes[0][0][0])* N_BLOCKS * ARRAY_SIZE * 3);
        for (h = 0; h < N_BLOCKS; h++)
        {
            gsl_vector_set_zero(IC_Stokes[0+h*3]);
            gsl_vector_set_zero(IC_Stokes[1+h*3]);
            gsl_vector_set_zero(IC_Stokes[2+h*3]);
        }
        memset(urad_array_perp, 0, sizeof(urad_array_perp[0][0]) * N_BLOCKS * ARRAY_SIZE); //reset these to 0 for next section
        memset(urad_array_para, 0, sizeof(urad_array_para[0][0]) * N_BLOCKS * ARRAY_SIZE);
        memset(dfactor_perp, 0, sizeof(dfactor_perp[0][0][0]) * N_BLOCKS * N_BLOCKS * ARRAY_SIZE); //reset these to 0 for next section
        memset(dfactor_para, 0, sizeof(dfactor_para[0][0][0]) * N_BLOCKS * N_BLOCKS * ARRAY_SIZE);

        //Opacity Stuff for integral after while loop
        dx_op_array[nSteps] = dx;

        nSteps += 1; //allows the number of jet sections to be determined
        PRINT("nStep: %d\n", nSteps);
        PRINT("x: %.5e\n", x);
        PRINT("B: %.5e\n", B);
        PRINT("R: %.5e\n", R);

        if (nSteps == breakstep)
        {
            break;
        }
    } // End of Jet calculation Loop ----------------------------------------------------------------------------------------//

    free(R_array);
    free(dx_array); //free all dynamically allocated memory
    for (i = 0; i < buffer_size; i++)
    {
        free(Pperp_array[i]);
        free(Ppara_array[i]);
    }
    free(Pperp_array);
    free(Ppara_array);

    for (h = 0; h < N_BLOCKS; h++)
    {
        gsl_vector_free(X[h]); //perp and para concatenated

        gsl_matrix_free(A[0+h*3]); //one for each stokes parameter
        gsl_matrix_free(A[1+h*3]);
        gsl_matrix_free(A[2+h*3]);

        gsl_vector_free(IC_Stokes[0+3*h]);
        gsl_vector_free(IC_Stokes[1+3*h]);
        gsl_vector_free(IC_Stokes[2+3*h]);
    }
    

    //Opacity integral calculation:
    for (n = 0; n < breakstep; n++)
    {
        for (l = 0; l < ARRAY_SIZE; l++)
        {
            //printf("K %d\t%d\t%.5e\n", n, l, k_array[n][l]);
            for (i = 0; i < breakstep; i++)
            {
                if (i > n)
                {
                    tau[n][l] += k_array[i][l] * dx_op_array[i] / cos(deg2rad(theta_obs));
                }
            }
        }
    }

    for (i = 0; i < breakstep; i++)
    {
        for (l = 0; l < ARRAY_SIZE; l++)
        {
            S_StokesTotal[l + 0 * ARRAY_SIZE] += S_StokesTotal_Sec[i][l][0] * exp(-tau[i][l]); // not quite right since binshifts, but approx ok for purposes of fit
            S_StokesTotal[l + 1 * ARRAY_SIZE] += S_StokesTotal_Sec[i][l][1] * exp(-tau[i][l]); //should address binshifts properly tho as it affects polarization of mm and below
            S_StokesTotal[l + 2 * ARRAY_SIZE] += S_StokesTotal_Sec[i][l][2] * exp(-tau[i][l]);
            //printf("TAU %d\t%d\t%.5e\n", i, l, tau[i][l]);
        }
    }

    for (n = 0; n < ARRAY_SIZE; n++)
    { //multiplying F(nu) by nu to get nuF(nu), to plot Power instead: switch f_pol by dfreqs_pol, however this will require a different treatment of ICS_total

        IC_StokesTotal[n + 0 * ARRAY_SIZE] *= f_pol_IC[n];
        IC_StokesTotal[n + 1 * ARRAY_SIZE] *= f_pol_IC[n];
        IC_StokesTotal[n + 2 * ARRAY_SIZE] *= f_pol_IC[n];

        S_StokesTotal[n + 0 * ARRAY_SIZE] *= f_pol[n];
        S_StokesTotal[n + 1 * ARRAY_SIZE] *= f_pol[n];
        S_StokesTotal[n + 2 * ARRAY_SIZE] *= f_pol[n];
    }

    return 1;
}

//executable test for debugging
int main(){
    double argv[] = {pow(10,38.436129), pow(10,10.12309571), 1.7230574,  39.97781039,
               15.28831413, pow(10,-4.39991182), 3.4409881, 1.07342862, pow(10,5.76709691)};
    int argblocks[] = {1, 0, 42, 0};
    double IC_StokesTotal[50] = {0};
    double S_StokesTotal[50] = {0};
    double f_pol_IC[50] = {0}; 
    double f_pol[50] = {0};

    int res = jetmodel(argv, argblocks, IC_StokesTotal, S_StokesTotal, f_pol_IC, f_pol);
    return res;
} 