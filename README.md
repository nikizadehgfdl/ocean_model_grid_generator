# GFDL Ocean Model Grid Generator
A collection of tools for creating finite element spherical tripolar grids for GFDL's MOM based ocean models.

Required python packages:
- numpypi : python -m pip install git+https://github.com/underwoo/numpypi@pip.installable
 
To test this module quickly try
- cd extras ; make -f Makefile.examples quick

Examples:
- To build a grid consisting a regular spherical grid between 50S and 70S, stitched with  a displaced pole cap south of 70S with the south pole displaced to (60E,50S) 
ocean_grid_generator.py -f 4 --south_cutoff_row 22 --write_subgrid_files --plot --gridlist=sc,so  --no-metrics --londp 60.0 --latdp -85.85702261518855 --south_cap_lat -70. --south_ocean_upper_lat -50.

[Technical guide](https://github.com/nikizadehgfdl/grid_generation/blob/dev/ocean_grid_generator_guide.pdf)

[![Build Status](https://travis-ci.org/nikizadehgfdl/ocean_model_grid_generator.svg?branch=dev)](https://travis-ci.org/nikizadehgfdl/ocean_model_grid_generator)

