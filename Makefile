TARGS = tripolar_res2.nc \
        tripolar_disp_res4.nc \
	tripolar_disp_res8.nc

all: $(TARGS) hash.md5
	md5sum -c hash.md5
       
tripolar_res2.nc: 
	./ocean_grid_generator.py -f tripolar_res2.nc -r 2 --rdp 0
tripolar_disp_res4.nc:
	./ocean_grid_generator.py -f tripolar_disp_res4.nc -r 4 --trim_south_80
tripolar_disp_res8.nc:
	./ocean_grid_generator.py -f tripolar_disp_res8.nc -r 8
       

hash.md5: | $(TARGS)
	md5sum $(TARGS) > $@

clean:
	rm -f $(TARGS) $(DEPS) pickle.*        
