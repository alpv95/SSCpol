# SSCpol
Python + C codebase to calculate the observed steady-state SED and polarization of a relativistic conical jet of electrons & positrons. Synchrotron and Synchrotron-Self Compton (SSC) Polarization are included.

Developed for [Peirson & Romani (2018)](https://iopscience.iop.org/article/10.3847/1538-4357/aad69d/meta), [Peirson & Romani (2019)](https://iopscience.iop.org/article/10.3847/1538-4357/ab46b1/meta) & [Peirson, Liodakis & Romani (2022)]; model inspired partly by  [Potter & Cotter (2012)](https://academic.oup.com/mnras/article/423/1/756/1747479) & [TEMZ](https://www.bu.edu/blazars/temz.html).

Model free parameters are: 
* jet power [W], 
* initial electron power law index, 
* the electron exponential cutoff energy [eV], 
* initial magnetic flux density at jet base [T], 
* bulk Lorentz factor, 
* observation angle [deg], 
* jet opening angle [deg], 
* initial equipartition fraction, 
* number of random B-field zones in one jet cross section (either 1,7,19,37,61,91,127). 

# Requirements
Python requirements noted in `requirements.txt`.
C requirements are noted in `Makefile`. (gcc, GSL, CBlas, see https://www.gnu.org/software/gsl/)

# Installation
git clone the repository and run `make` in `src/` to build the C code.

# Usage
See `example.ipynb` for common usage.
`python3 -u fit.py` to run the model fits to Blazar spectra.
`plot_results.ipynb` to visualize results and reproduce plots from [Peirson, Liodakis & Romani (2022)].
`sscpol/jet_fns.py: run_ssc()` for single model evaluation.

# Attribution
```
    @article{peirson_polarization_2019,
	title = {The {Polarization} {Behavior} of {Relativistic} {Synchrotron} {Self}-{Compton} {Jets}},
	volume = {885},
	issn = {0004-637X},
	url = {https://doi.org/10.3847%2F1538-4357%2Fab46b1},
	doi = {10.3847/1538-4357/ab46b1},
	language = {en},
	number = {1},
	urldate = {2020-03-28},
	journal = {ApJ},
	author = {Peirson, A. L. and Romani, Roger W.},
	month = nov,
	year = {2019},
	pages = {76},
}
```

Basic algorithm:
-----
Numbers in paretheses refer to equations in [Peirson & Romani (2019)].
![Flow chart outlining the basic algorithm](FlowChart.jpg?raw=true "Title")