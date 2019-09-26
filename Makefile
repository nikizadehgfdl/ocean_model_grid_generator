TARGS = ocean_hgrid_res4.0.nc \
        ocean_hgrid_res1.0.nc \
        ocean_hgrid_res0.5.nc \
        ocean_hgrid_res0.5_equenh.nc \
        ocean_hgrid_res0.25.nc \
        ocean_hgrid_res0.125.nc 

#Note: Github Travis cannot make the higher res grids below and errors out with "MemoryError"
#      That is why they are commented out from TARGS so Travis can finish to keep records.   
QUICKTARGS = ocean_hgrid_res4.0.nc \
        ocean_hgrid_res1.0.nc \
        ocean_hgrid_res0.5.nc \
        ocean_hgrid_res0.5_equenh.nc

all: $(TARGS) hash.md5
	cat hash.md5
	md5sum -c hash.md5

quick : $(QUICKTARGS) 
	head -4 hash.md5 > hash.quick
	md5sum -c hash.quick

#Note: --no_changing_meta arg is used to avoid putting time/platform dependent info in the files so that they can  be checksumed.
#      Please do not use this arg for normal grid generation, it prevents adding useful information to meta data.
ocean_hgrid_res4.0.nc:
	time ./ocean_grid_generator.py $(DEBUG) -f ocean_hgrid_res4.0.nc -r 0.25 --rdp 0 --no_changing_meta
ocean_hgrid_res1.0.nc:
	time ./ocean_grid_generator.py $(DEBUG) -f ocean_hgrid_res1.0.nc -r 1.0  --rdp 0 --south_cutoff_row 2 --no_changing_meta
ocean_hgrid_res0.5.nc: 
	time ./ocean_grid_generator.py $(DEBUG) -f ocean_hgrid_res0.5.nc -r 2    --rdp 0 --no_changing_meta 
ocean_hgrid_res0.5_equenh.nc: 
	time ./ocean_grid_generator.py $(DEBUG) -f ocean_hgrid_res0.5_equenh.nc -r 2 --rdp 0 --south_cutoff_row 128 --no_changing_meta --write_subgrid_files --enhanced_equatorial
ocean_hgrid_res0.25.nc:
	time ./ocean_grid_generator.py $(DEBUG) -f ocean_hgrid_res0.25.nc -r 4 --south_cutoff_row 83 --write_subgrid_files --no_changing_meta
ocean_hgrid_res0.125.nc:
	time ./ocean_grid_generator.py $(DEBUG) -f ocean_hgrid_res0.125.nc -r 8 --south_cutoff_row 5 --write_subgrid_files --no_changing_meta


hash.md5: | $(TARGS)
	md5sum $(TARGS) > $@
	cat $@

check:
	md5sum -c hash.md5

clean:
	rm -f $(TARGS) $(DEPS) ocean_hgrid_res*.nc
