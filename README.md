[![DOI](https://zenodo.org/badge/105955957.svg)](https://zenodo.org/badge/latestdoi/105955957)

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
* Python requirements noted in `requirements.txt` and automatically installed when running installation.
* C requirements are noted in `Makefile`. (gcc, GSL, CBlas, see https://www.gnu.org/software/gsl/)

# Installation
git clone this repository then `pip3 install -e .` from root directory.

# Usage
* See `notebooks/example.ipynb` for common usage.
* `python3 -u fit.py` to run the model fits to Blazar spectra. 
* See `python3 fit.py -h`
* `notebooks/plot_results.ipynb` to visualize results and reproduce plots from [Peirson, Liodakis & Romani (2022)].
`from sscpol.jet_fns import run_ssc` for single model evaluation.

# Testing
Once installed run: `python3 -m pytest sscpol` to execute the test suite.

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
