# GFDL Ocean Model Grid Generator
A collection of tools for creating finite element spherical tripolar grids for GFDL's MOM based ocean models.

Required python packages:
- numpypi : python -m pip install git+https://github.com/underwoo/numpypi@pip.installable
 
To test this module quickly try
- export PYTHONPATH="/net2/nnz/opt/miniconda/lib/python3.8/site-packages"; export PATH="/net2/nnz/opt/miniconda/bin:$PATH"; export HDF5_USE_FILE_LOCKING=FALSE ; . activate om4labs
- cd extras ; make -f Makefile.examples quick

[Technical guide](https://github.com/nikizadehgfdl/grid_generation/blob/dev/ocean_grid_generator_guide.pdf)

[![Build Status](https://travis-ci.org/nikizadehgfdl/ocean_model_grid_generator.svg?branch=dev)](https://travis-ci.org/nikizadehgfdl/ocean_model_grid_generator)

