TARGS = tripolar_res2.nc \
        tripolar_res2_MIDAS.nc \
        tripolar_disp_res4.nc \
        tripolar_disp_res4_MIDAS.nc \
	tripolar_res8_repold8.nc \
	tripolar_disp_res8.nc

all: $(TARGS) hash.md5
	cat hash.md5
	md5sum -c hash.md5
       
#Note: --no_changing_meta arg is used to avoid putting time/platform dependent info in the files so that they can  be checksumed.
#      Please do not use this arg for normal grid generation, it prevents adding useful information to meta data.
tripolar_res2.nc: 
	./ocean_grid_generator.py -f tripolar_res2.nc -r 2 --rdp 0 --no_changing_meta 
tripolar_res2_MIDAS.nc:
	./ocean_grid_generator.py -f tripolar_res2_MIDAS.nc -r 2 --rdp 0 --south_cutoff_row 81 --reproduce_MIDAS_grids --no_changing_meta
tripolar_disp_res4_MIDAS.nc:
	./ocean_grid_generator.py -f tripolar_disp_res4_MIDAS.nc -r 4 --south_cutoff_row 81 --reproduce_MIDAS_grids --no_changing_meta
tripolar_disp_res4.nc:
	./ocean_grid_generator.py -f tripolar_disp_res4.nc -r 4 --south_cutoff_row 81 --no_changing_meta
tripolar_disp_res8.nc:
	./ocean_grid_generator.py -f tripolar_disp_res8.nc -r 8 --write_subgrid_files --no_changing_meta
tripolar_res8_repold8.nc:
	./ocean_grid_generator.py -f tripolar_res8_repold8.nc -r 8 --reproduce_old8_grids --no_changing_meta
       

hash.md5: | $(TARGS)
	md5sum $(TARGS) > $@
	cat $@

check:
	md5sum -c hash.md5

clean:
	rm -f $(TARGS) $(DEPS) pickle.*        
