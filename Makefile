TARGS = tripolar_res2.nc \
        tripolar_res2_MIDAS.nc \
        tripolar_disp_res4.nc \
        tripolar_disp_res4_MIDAS.nc \
	tripolar_res8_repold8.nc \
	tripolar_disp_res8.nc

all: $(TARGS) hash.md5.githubTravisCI
	md5sum -c hash.md5.githubTravisCI 
       
tripolar_res2.nc: 
	./ocean_grid_generator.py -f tripolar_res2.nc -r 2 --rdp 0 
tripolar_res2_MIDAS.nc:
	./ocean_grid_generator.py -f tripolar_res2_MIDAS.nc -r 2 --rdp 0 --south_cutoff_row 81 --reproduce_MIDAS_grids
tripolar_disp_res4_MIDAS.nc:
	./ocean_grid_generator.py -f tripolar_disp_res4_MIDAS.nc -r 4 --south_cutoff_row 81 --reproduce_MIDAS_grids
tripolar_disp_res4.nc:
	./ocean_grid_generator.py -f tripolar_disp_res4.nc -r 4 --south_cutoff_row 81 
tripolar_disp_res8.nc:
	./ocean_grid_generator.py -f tripolar_disp_res8.nc -r 8 --write_subgrid_files
tripolar_res8_repold8.nc:
	./ocean_grid_generator.py -f tripolar_res8_repold8.nc -r 8 --reproduce_old8_grids
       

hash.md5.githubTravisCI: | $(TARGS)
	md5sum $(TARGS) > $@
	cat $@

clean:
	rm -f $(TARGS) $(DEPS) pickle.*        
