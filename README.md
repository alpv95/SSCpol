# Synchrotron and SSC Polarization of Blazars
Python and C code to calculate the SED and polarization of relativistic conical jet of electrons and positrons initially in equipartition with power law $\alpha$ and an exponential cutoff. Model described in apj1 and apj2. The free parameters are: jet power  $W_j$, initial electron power law $\alpha$, the electron energy exponential cutoff $E_{max}$,  initial magnetic flux density $B_0$, relativisic bulk gamma $\Gamma_{bulk}$, observation angle $\theta_{obs}$, jet opening angle (jet frame) $\theta_{op}$ and the number of random B-field zones in one jet cross section $n_blocks$ (1,7,19,37,61,91,127). Isotropic random B-fields are generated using the Mersenne twister algorithm.

## To run
jet_model.c contains the full simulation, outputting a number of text files with the results: basicdata.txt (jet dx, x, B and R at each simulation step), keyparams.txt (jet initial parameters), pi.txt (polarization fraction, EVPA and observed power results for synchrotron, SSC and joint) and IC_Z.txt / S_Z.txt (Stokes parameters for individual jet B-field zones).

runplotSED.py acts as a wrapper for jet_model.c and allows one to run and plot the results. Use --help to see options and flags. Jet parameters other than $\theta_{obs}$ can be modified in jet_model.c directly.

![Alt text](SEDsingle.pdf?raw=true "Title")





