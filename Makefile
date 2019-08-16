TARGS = ocean_hgrid_res4.0.nc \
        ocean_hgrid_res1.0.nc \
        ocean_hgrid_res0.5.nc \
        ocean_hgrid_res0.5_equenh.nc \
        ocean_hgrid_res0.25.nc \
        ocean_hgrid_res0.125.nc 

MIDAG = ocean_hgrid_res0.5_MIDAS.nc \
	ocean_hgrid_res0.25_MIDAS.nc 

MIDAS: $(MIDAG) check

all: $(TARGS) hash.md5
	cat hash.md5
	md5sum -c hash.md5

#Note: --no_changing_meta arg is used to avoid putting time/platform dependent info in the files so that they can  be checksumed.
#      Please do not use this arg for normal grid generation, it prevents adding useful information to meta data.
ocean_hgrid_res4.0.nc:
	time ./ocean_grid_generator.py $(DEBUG) -f ocean_hgrid_res4.0.nc -r 0.25 --rdp 0 --no_changing_meta
ocean_hgrid_res1.0.nc:
	time ./ocean_grid_generator.py $(DEBUG) -f ocean_hgrid_res1.0.nc -r 1.0  --rdp 0 --no_changing_meta
ocean_hgrid_res0.5.nc: 
	time ./ocean_grid_generator.py $(DEBUG) -f ocean_hgrid_res0.5.nc -r 2    --rdp 0 --no_changing_meta 
ocean_hgrid_res0.25.nc:
	time ./ocean_grid_generator.py $(DEBUG) -f ocean_hgrid_res0.25.nc -r 4 --south_cutoff_row 83 --write_subgrid_files --no_changing_meta
ocean_hgrid_res0.125.nc:
	time ./ocean_grid_generator.py $(DEBUG) -f ocean_hgrid_res0.125.nc -r 8 --south_cutoff_row 5 --write_subgrid_files --no_changing_meta
	#module swap python/3.6.4 On gfdl pan
ocean_hgrid_res0.5_equenh.nc: 
	time ./ocean_grid_generator.py $(DEBUG) -f ocean_hgrid_res0.5_equenh.nc -r 2 --rdp 0 --south_cutoff_row 128 --no_changing_meta --write_subgrid_files --enhanced_equatorial
ocean_hgrid_res0.25_smooth.nc:
	time ./ocean_grid_generator.py $(DEBUG) -f ocean_hgrid_res0.25_smooth.nc -r 4 --south_cutoff_row 65 --smooth_dy --write_subgrid_files --no_changing_meta
ocean_hgrid_res0.125_smooth.nc:
	time ./ocean_grid_generator.py $(DEBUG) -f ocean_hgrid_res0.125_smooth.nc -r 8 --south_cutoff_row 5 --smooth_dy --write_subgrid_files --no_changing_meta
	#module swap python/3.6.4 On gfdl pan
ocean_hgrid_res0.5_equenh_smooth.nc: 
	time ./ocean_grid_generator.py $(DEBUG) -f ocean_hgrid_res0.5_equenh_smooth.nc -r 2 --rdp 0 --south_cutoff_row 122 --smooth_dy --no_changing_meta --write_subgrid_files --enhanced_equatorial


#MIDAS grids, work only on gfdl PAN
ocean_hgrid_res0.5_MIDAS.nc:
	#module swap python/2.7.3_workstation On gfdl pan
	time ./ocean_grid_generator.py $(DEBUG) -f ocean_hgrid_res0.5_MIDAS.nc -r 2 --rdp 0 --south_cutoff_row 130 --reproduce_MIDAS_grids --write_subgrid_files --no_changing_meta
	#module load nccmp 
	nccmp -d /archive/gold/datasets/OM4_05/mosaic_ocean.v20180227.unpacked/ocean_hgrid.nc ocean_hgrid_res0.5_MIDAS.nc
ocean_hgrid_res0.25_MIDAS.nc:
	#module swap python/2.7.3_workstation On gfdl pan
	time ./ocean_grid_generator.py $(DEBUG) -f ocean_hgrid_res0.25_MIDAS.nc -r 4 --south_cutoff_row 81 --reproduce_MIDAS_grids --write_subgrid_files --no_changing_meta
	#module load nccmp 
	nccmp -d /archive/gold/datasets/OM4_025/mosaic.v20170622.unpacked/ocean_hgrid.nc ocean_hgrid_res0.25_MIDAS.nc

hash.md5: | $(TARGS)
	md5sum $(TARGS) > $@
	cat $@

check:
	md5sum -c hash.md5

clean:
	rm -f $(TARGS) $(DEPS) ocean_hgrid_res*.nc
