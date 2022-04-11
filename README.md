# Synchrotron and Synchrotron-Self Compton (SSC) Polarization of Blazars
Codebase to calculate the observed SED and polarization of relativistic conical jet of electrons and positrons initially in equipartition following a power law with an exponential cutoff energy. Developed for [Peirson & Romani (2018)](https://iopscience.iop.org/article/10.3847/1538-4357/aad69d/meta) & [Peirson & Romani (2019)]; model inspired partly by  [Potter & Cotter (2012)](https://academic.oup.com/mnras/article/423/1/756/1747479) & [TEMZ](https://www.bu.edu/blazars/temz.html).
The free parameters are: jet power [W], initial electron power law index, the electron exponential cutoff energy [J],  initial magnetic flux density at jet base [T], bulk Lorentz factor, observation angle [rad], jet opening angle [rad] and the number of random B-field zones in one jet cross section (either 1,7,19,37,61,91,127). Isotropic random B-fields are generated using the Mersenne twister algorithm (cite). The SSC calculation is parallelized using **OpenMP**.
# SSC Fast
Sped up version of SSC_final using the magic of linear algebra. Requires GSL C library (can be download and installation instructions here: https://www.gnu.org/software/gsl/). 

##Shared library (.so) build
gcc -fPIC -shared func_model.c mtwister.c mtwister.h jet_fns.c jet_fns.h -lgsl -lgslblas -lm

# Synchrotron and SSC Polarization of Blazars
Python and C code to calculate the SED and polarization of relativistic conical jet of electrons and positrons initially in equipartition with power law $\alpha$ and an exponential cutoff. Model described in apj1 and apj2. The free parameters are: jet power  $W_j$, initial electron power law $\alpha$, the electron energy exponential cutoff $E_{max}$,  initial magnetic flux density $B_0$, relativisic bulk gamma $\Gamma_{bulk}$, observation angle $\theta_{obs}$, jet opening angle (jet frame) $\theta_{op}$ and the number of random B-field zones in one jet cross section $n_blocks$ (1,7,19,37,61,91,127). Isotropic random B-fields are generated using the Mersenne twister algorithm.

## To run
jet_model.c contains the full simulation, outputting a number of text files with the results: `basicdata.txt` (jet *dx*, *x*, *B* and *R* at each simulation step), `keyparams.txt` (jet initial parameters), `pi.txt` (polarization fraction, EVPA and observed power results for synchrotron, SSC and joint) and `IC_Z.txt` / `S_Z.txt` (Stokes parameters for individual jet B-field zones).

`runplotSED.py` acts as a wrapper for jet_model.c and allows one to run and plot the results. Use `--help` to see options and flags. For example to run a full sync+SSC jet model with 7 B-field blocks for 100 steps at an observation angle of 1.5 deg:
``` 
python runplotSED.py 1.5 -ssc --nblocks 7 --nsteps 100
```
Jet parameters other than the observation angle can be modified in `jet_model.c` directly. The folder `utils` contains a number of plotting and data analysis functions to display the results of the model.

## Requirements
Python 3 with and gcc are the only requirements (gcc version requires OpenMP support if parallelization is desired).


Basic algorithm:
-----
Numbers in paretheses refer to equations in [Peirson & Romani (2019)].
![Flow chart outlining the basic algorithm](FlowChart.jpg?raw=true "Title")




