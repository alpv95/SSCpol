# from setuptools import find_packages, setup

# setup(
#     name="sscpol",
#     packages=find_packages(where="."),
#     package_dir={"": "."},
# )

from setuptools import setup
from distutils.core import Extension
ext_fast = Extension(name='sscpol.func_model', 
                    include_dirs = ['include'],
                    libraries = ['m', 'gsl', 'gslcblas'],
                    # library_dirs = ['/usr/local/lib'],
                    sources=['src/func_model.c','src/jet_fns.c','src/mtwister.c'])
# ext_full = Extension(name='sscpol.full_model', 
#                     include_dirs = ['include'],
#                     libraries = ['m', 'gsl', 'gslcblas'],
#                     sources=['src/full_model.c','src/jet_fns.c','src/mtwister.c'])
setup(use_scm_version=True, ext_modules=[ext_fast])