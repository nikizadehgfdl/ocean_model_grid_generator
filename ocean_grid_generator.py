#!/usr/bin/env python

import numpy as np
import sys, getopt
import datetime, os, subprocess

#Constants
PI_180 = np.pi/180.
#_default_Re = 6.378e6
_default_Re = 6371.e3 #MIDAS

def generate_bipolar_cap_grid(Ni,Nj_ncap,lat0_bp,lon_bp,lenlon):
    print( 'Generating bipolar grid bounded at latitude ',lat0_bp  )
    rp=np.tan(0.5*(90-lat0_bp)*PI_180)
    #First define a (lon,lat) coordinate on the Northern hemisphere of the globe sphere
    #such that the resolution of latg matches the desired resolution of the final grid along the symmetry meridian 
    lon_g = lon_bp  + np.arange(Ni+1) * lenlon/Ni 
    lamg = np.tile(lon_g,(Nj_ncap+1,1)) 
    latg0_cap = lat0_bp + np.arange(Nj_ncap+1) * (90-lat0_bp)/Nj_ncap
    phig0_cap = np.tile(latg0_cap.reshape((Nj_ncap+1,1)),(1,Ni+1))
    ### symmetry meridian resolution fix 
    phig = 90-2*np.arctan(np.tan(0.5*(90-phig0_cap)*PI_180)/rp)/PI_180

    #Simplify  the formulas to avoid division by zero
    #alpha  = np.cos((lamg-lon_p)*PI_180) 
    alpha2 = (np.cos((lamg-lon_bp)*PI_180))**2
    #beta = -np.cotan(phig*PI_180)
    beta2_inv = (np.tan(phig*PI_180))**2
    
    A=np.sqrt(1-alpha2)*np.sin(phig*PI_180) #Actually two equations  +- |A|    
    B=np.sqrt((1-alpha2)/(1+alpha2*beta2_inv)) #Actually two equations  +- |B|
#   Equivalently we can do the following which has manifest symmetry lam --> 180+lam
#    A=np.sin((lamg-lon_bp)*PI_180)*np.sin(phig*PI_180) #Actually two equations  +- |A|
#    A=np.where((lamg-lon_bp)>180,-A,A)
#    
#    B=np.sin((lamg-lon_bp)*PI_180)/np.sqrt(1+alpha2*beta2_inv) #Actually two equations  +- |B|
#    B=np.where((lamg-lon_bp)>180,-B,B)

    #Deal with beta=0
    B=np.where(np.abs(beta2_inv)>1.0E10 , 0.0, B)

    lamc = np.arcsin(B)/PI_180 
    #phic = np.arcsin(A)/PI_180
    #or use
    chic = np.arccos(A)

    ##But this equation accepts 4 solutions for a given B, {l, 180-l, l+180, 360-l } 
    ##We have to pickup the "correct" root. 
    ##One way is simply to demand lamc to be continuous with lam on the equator phi=0
    ##I am sure there is a more mathematically concrete way to do this.
    lamc = np.where((lamg-lon_bp>90)&(lamg-lon_bp<=180),180-lamc,lamc)
    lamc = np.where((lamg-lon_bp>180)&(lamg-lon_bp<=270),180+lamc,lamc)
    lamc = np.where((lamg-lon_bp>270),360-lamc,lamc)
    #Along symmetry meridian choose lamc
    lamc = np.where((lamg-lon_bp==90),90,lamc)    #Along symmetry meridian choose lamc=90-lon_bp
    lamc = np.where((lamg-lon_bp==270),270,lamc)  #Along symmetry meridian choose lamc=270-lon_bp    
    lams = lamc + lon_bp

    ##Project back onto the larger (true) sphere so that the projected equator shrinks to latitude \phi_P=lat0_tp
    ##then we have tan(\phi_s'/2)=tan(\phi_p'/2)tan(\phi_c'/2)

    #phis = 90 - 2 * np.arctan(rp * np.tan(0.5*(90-phic)*PI_180))/PI_180
    #or equivalently
    phis = 90 - 2 * np.arctan(rp * np.tan(chic/2))/PI_180
    print('   number of js=',phis.shape[0])

    ##Calculate the Metrics
    M_inv = rp * (1 + (np.tan(chic/2))**2) / (1 + (rp*np.tan(chic/2))**2)
    chig = (90-phig)*PI_180
    N     = rp * (1 + (np.tan(chig/2))**2) / (1 + (rp*np.tan(chig/2))**2)
    N_inv = 1/N    
    cos2phis = (np.cos(phis*PI_180))**2 

    h_j_inv = cos2phis*alpha2*(1-alpha2)*beta2_inv*(1+beta2_inv)/(1+alpha2*beta2_inv)**2 \
            +  M_inv*M_inv*(1-alpha2)/(1+alpha2*beta2_inv) 
    h_j_inv = np.sqrt(h_j_inv)*N_inv*(90-lat0_bp)*PI_180/Nj_ncap
#    h_j_inv = h_j_inv[:-1,:]

    h_i_inv = cos2phis * (1+beta2_inv)/(1+alpha2*beta2_inv)**2 + M_inv*M_inv*alpha2*beta2_inv/(1+alpha2*beta2_inv)
    h_i_inv = np.sqrt(h_i_inv)*(2*np.pi/Ni)
#    h_i_inv = h_i_inv[:,:-1]

    return lams,phis,h_i_inv,h_j_inv

def bp_lam(x,y,bpeq,rp):
    """bp_lam = ((90-y)/(90-lat_join))*90
       invert Murray's eqn. 5b with phic=0 to place point at specified geo. lat """
    bp_lam = 2.*np.arctan(np.tan((0.5*np.pi-y*PI_180)/2)/rp)/PI_180
    bp_lam = np.where(mdist(x,bpeq)<90.,-bp_lam, bp_lam)
    return bp_lam    

def bp_phi(x,y,bpsp,bpnp):
    bps = mdist(x,bpsp)
    bpn = mdist(x,bpnp)
    bp_phi = np.where(bps<90,-90+bps,90-bpn)
    return bp_phi

def generate_bipolar_cap_grid_fms(Ni,Nj_ncap,lat0_p,lon_p,lenlon,lenlat):
    print( 'Generating FMS bipolar grid bounded at latitude ',lat0_p  )
    rp=np.tan(0.5*(90-lat0_p)*PI_180) #Murray section 2.2 before Eq(6) r_p=tan(\phi_P\prime /2) 
                                                   #where \phi_P is the latitude of the bounding parrallel lat0
    sp_lon_cap = lon_p  + np.arange(Ni+1) * lenlon/Ni 
    sp_lat_cap = lat0_p + np.arange(Nj_ncap+1) * lenlat/Nj_ncap
    lon_cap = np.tile(sp_lon_cap,(Nj_ncap+1,1)) 
    lat_cap = np.tile(sp_lat_cap.reshape((Nj_ncap+1,1)),(1,Ni+1))

    lamc_fms = bp_lam(lon_cap,lat_cap,lon_p+90,rp) 
    phic_fms = bp_phi(lon_cap,lat_cap,lon_p,lon_p+180)
    lams_fms = lon_p + 90*0 - np.arctan2(np.sin(lamc_fms*PI_180),np.tan(phic_fms*PI_180))/PI_180 #eqn.5a
    #Note the *0:          |
    #The original grid is rotated 90 degrees compared to both MIDAS and new

    chic_fms = np.arccos(np.cos(lamc_fms*PI_180)*np.cos(phic_fms*PI_180)) #eqn.6
    phis_fms = 90 - 2 * np.arctan(rp*np.tan(chic_fms/2))/PI_180
    print('   number of js=',phis_fms.shape[0])
    return lams_fms, phis_fms



def y_mercator(Ni, phi):
    """Equation (1)"""
    R = Ni / (2 * np.pi)
    return R * ( np.log( (1.0 + np.sin(phi) ) / np.cos(phi)) )
def phi_mercator(Ni, y):
    """Equation (2)"""
    R = Ni / (2 * np.pi)
    return np.arctan( np.sinh(y/R) ) * (180/np.pi) # Converted to degrees
def y_mercator_rounded(Ni, phi):
    y_float = y_mercator(Ni, phi)
    return ( np.sign(y_float) * np.ceil( np.abs(y_float) ) ).astype(int)

def generate_mercator_grid(Ni,phi_s,phi_n,lon0_M,lenlon_M,shift_equator_to_u_point=True, ensure_nj_even=True):
    print( 'Requesting Mercator grid with phi range: phi_s,phi_n=', phi_s,phi_n )
    # Diagnose nearest integer y(phi range)
    y_star = y_mercator_rounded(Ni, np.array([phi_s*PI_180,phi_n*PI_180]))
    print( '   y*=',y_star, 'nj=', y_star[1]-y_star[0]+1  )
    #Ensure that the equator (y=0) is a u-point
    if(y_star[0]%2 == 0):
        print("   Equator is not going to be a u-point!")
        if(shift_equator_to_u_point):
            print("   Fixing this by shifting the bounds!")
            y_star[0] = y_star[0] - 1
            y_star[1] = y_star[1] - 1
            print( '   y*=',y_star, 'nj=', y_star[1]-y_star[0]+1 )
    if((y_star[1]-y_star[0]+1)%2 == 0 and ensure_nj_even):
        print("   Supergrid has an odd number of area cells!")
        if(ensure_nj_even):
            print("   Fixing this by shifting the y_star[1] ")
            y_star[1] = y_star[1] - 1
            print( '   y*=',y_star, 'nj=', y_star[1]-y_star[0]+1 )
    Nj=y_star[1]-y_star[0]
    print( '   Generating Mercator grid with phi range: phi_s,phi_n=', phi_mercator(Ni, y_star) )
    phi_M = phi_mercator(Ni, np.arange(y_star[0],y_star[1]+1)) 
    #Ensure that the equator (y=0) is included and is a u-point
    equator=0.0
    equator_index = np.searchsorted(phi_M,equator)
    if(equator_index == 0): 
        raise Exception('   Ooops: Equator is not in the grid')
    else:
        print("   Equator is at j=", equator_index)
    #Ensure that the equator (y=0) is a u-point
    if(equator_index%2 == 0):
        raise Exception("Ooops: Equator is not going to be a u-point")

    y_grid_M = np.tile(phi_M.reshape(Nj+1,1),(1,Ni+1))
    lam_M = lon0_M + np.arange(Ni+1) * lenlon_M/Ni
    x_grid_M = np.tile(lam_M,(Nj+1,1)) 
    print('   number of js=',y_grid_M.shape[0])
    return x_grid_M,y_grid_M

                                
def generate_displaced_pole_grid(Ni,Nj_scap,lon0,lenlon,lon_dp,r_dp,lat0_SO,doughnut,nparts=8, ensure_nj_even=True):
    print( 'Generating displaced pole grid bounded at latitude ',lat0_SO  )
    print('   rdp=',r_dp,' , doughnut=',doughnut)
    x=lon0 + np.arange(Ni+1) * lenlon/Ni
    y=np.linspace(-90.,0.5*(lat0_SO-90.0),Nj_scap//nparts)
    y=np.concatenate((y,np.linspace(y.max(),lat0_SO,1+Nj_scap*(nparts-1)//nparts)))
    if(y.shape[0]%2 == 0 and ensure_nj_even):
        print("   The number of j's is not even. Fixing this by cutting one row at south.")
        y = np.delete(y,0,0)
    X1,Y1=np.meshgrid(x,y)
    lamc_DP,phic_DP = displaced_pole_cap(X1,Y1,lam_pole=-lon_dp,r_pole=r_dp,lat_joint=lat0_SO,
                                         excluded_fraction=doughnut) 
    print('   number of js=',phic_DP.shape[0])
    return lamc_DP,phic_DP


def displaced_pole_cap(lon_grid,lat_grid,lam_pole,r_pole,lat_joint,excluded_fraction=None):

    #Projection from center of globe to plane tangent at south pole
    r_joint = np.tan((90+lat_joint)*PI_180) 
    z_0= r_pole * np.exp(1j*lam_pole*PI_180) 

    r = np.tan((90+lat_grid) *PI_180)/r_joint

    #Find the theta that has matching resolution at the unit circle with longitude at the joint
    #This is a conformal transformation the unit circle (inverse to the one below)
    e2itheta = np.exp(1j*lon_grid*PI_180) 
    e2ithetaprime = (e2itheta + z_0)/(1. + np.conj(z_0)*e2itheta)

    #Conformal map to displace pole from r=0 to r=r_dispole
    z=r*e2ithetaprime
    w=(z-z_0)/(1-np.conj(z_0)*z)
    
    #Inverse projection from tangent plane back to sphere
    lamcDP = np.angle(w, deg=True)
    #np.angle returns a value in the interval (-180,180)
    #However the input grid longitude is in (-lon0,-lon0+360), e.g., (-300,60)
    #We should shift the angle to be in that interval (note lon_grid[0,-1] = 60)
    #lamcDP = np.where(lamcDP>lon_grid[0,-1],lamcDP-360,lamcDP)
    #This will fix the discontinuity at lon_grid<0 but makes lamcDP discontinous at the end of the interval near 60
    #This is because there are a number of points with lamcDP>60 just below lon_grid=60
    #We need a second condition to leave these points intact, one that works is:
    lamcDP = np.where((lamcDP>lon_grid[0,-1])&(lon_grid<0) ,lamcDP-360,lamcDP)
    #Niki: The second condition above is ad hoc. Work on a more elaborate condition to get rid of the  discontinuity at lon_grid>60
    #Another adhoc fix for single point correction
    lamcDP[-1,0] = lamcDP[-1,0]-360

    rw=np.absolute(w)
    phicDP = -90+np.arctan(rw*r_joint)/PI_180
    if excluded_fraction is not None:
        ny,nx = lon_grid.shape 
        jmin=np.ceil(excluded_fraction*ny)
        jmin=jmin+np.mod(jmin,2)
        jmint = int(jmin)
        return lamcDP[jmint:,:], phicDP[jmint:,:]
    else:
        return lamcDP,phicDP

def cut_below(lam,phi,lowerlat):
    nj,ni = lam.shape
    for j in range(0,nj):
        if(phi[j,0]>lowerlat):
            break
    jmin=j
#    print("jmin",jmin)
    return lam[jmin:,:], phi[jmin:,:]

def cut_above(lam,phi,upperlat):
    nj,ni = lam.shape
    for j in range(0,nj):
        if(phi[j,0]>upperlat):
            break
    jmax=j
#    print("jmax",jmax)
    return lam[0:jmax,:], phi[0:jmax,:]

#utility function to plot grids
def plot_mesh_in_latlon(lam, phi, stride=1, phi_color='k', lam_color='r', newfig=True, title=None):
    import matplotlib.pyplot as plt
#    import seaborn as sns; sns.set()
    if (phi.shape != lam.shape): raise Exception('Ooops: lam and phi should have same shape')
    nj,ni = lam.shape
    if(newfig):
        plt.figure(figsize=(10,10))
    for i in range(0,ni,stride):
        plt.plot(lam[:,i],phi[:,i],lam_color)
    for j in range(0,nj,stride):
        plt.plot(lam[j,:],phi[j,:],phi_color)
    if title is not None:
        plt.title(title)
    plt.show()
        
def plot_mesh_in_xyz(lam, phi, stride=1, phi_color='k', lam_color='r', lowerlat=None, upperlat=None, newfig=True, title=None):
    if lowerlat is not None:
        lam,phi = cut_below(lam,phi,lowerlat=lowerlat)        
    if upperlat is not None:
        lam,phi = cut_above(lam,phi,upperlat=upperlat)        
    x = np.cos(phi*PI_180) * np.cos(lam*PI_180)
    y = np.cos(phi*PI_180) * np.sin(lam*PI_180)
    z = np.sin(phi*PI_180)        
    plot_mesh_in_latlon(x, y, stride=stride, phi_color=phi_color, lam_color=lam_color, newfig=newfig, title=title)

def mdist(x1,x2):
  """Returns positive distance modulo 360."""
  a=np.mod(x1-x2+720.,360.)
  b=np.mod(x2-x1+720.,360.)
  d=np.minimum(a,b)
  return d

def generate_grid_metrics(x,y,axis_units='degrees',Re=_default_Re, latlon_areafix=False):
    nytot,nxtot = x.shape
    if  axis_units == 'm':
      metric=1.0
    if  axis_units == 'km':            
      metric=1.e3
    if  axis_units == 'degrees':                        
      metric=Re*PI_180
    ymid_j = 0.5*(y+np.roll(y,shift=-1,axis=0))
    ymid_i = 0.5*(y+np.roll(y,shift=-1,axis=1))      
    dy_j = np.roll(y,shift=-1,axis=0) - y
    dy_i = np.roll(y,shift=-1,axis=1) - y
    dx_i = mdist(np.roll(x,shift=-1,axis=1),x)
    dx_j = mdist(np.roll(x,shift=-1,axis=0),x)
    dx = metric*metric*(dy_i*dy_i + dx_i*dx_i*np.cos(ymid_i*PI_180)*np.cos(ymid_i*PI_180))
    dx = np.sqrt(dx)
    dy = metric*metric*(dy_j*dy_j + dx_j*dx_j*np.cos(ymid_j*PI_180)*np.cos(ymid_j*PI_180))
    dy = np.sqrt(dy)
    dx=dx[:,:-1]
    dy=dy[:-1,:]
    if(latlon_areafix):
        delsin_j = np.roll(np.sin(y*PI_180),shift=-1,axis=0) - np.sin(y*PI_180)
        area=metric*metric*dx_i[:-1,:-1]*delsin_j[:-1,:-1]/PI_180
    else:
        area=dx[:-1,:]*dy[:,:-1]    
    angle_dx=np.zeros((nytot,nxtot))
#    angle_dx = np.arctan2(dy_i,dx_i)/PI_180      
#    self.angle_dx = numpy.arctan2(dy_i,dx_i)*180.0/numpy.pi
    # The commented out code above was incorrect for non-Cartesian grids
    # The corrected version, in addition to including spherical metrics, is centered in the interior and one-sided at the grid edges
    angle_dx[:,1:-1] = np.arctan2(y[:,2:]-y[:,:-2],(x[:,2:]-x[:,:-2])*np.cos(y[:,1:-1]*PI_180))
    angle_dx[:,0]    = np.arctan2(y[:,1] -y[:,0]  ,(x[:,1] -x[:,0]  )*np.cos(y[:,0]*PI_180))
    angle_dx[:,-1]   = np.arctan2(y[:,-1]-y[:,-2] ,(x[:,-1]-x[:,-2] )*np.cos(y[:,-1]*PI_180))
    angle_dx = angle_dx /PI_180
    return dx,dy,area,angle_dx

def metrics_error(dx_,dy_,area_,Ni,lat1,lat2=90,Re=_default_Re):
    exact_area = 2*np.pi*(Re**2)*np.abs(np.sin(lat2*PI_180)-np.sin(lat1*PI_180))
    exact_lat_arc_length = np.abs(lat2-lat1)*PI_180*Re  
    exact_lon_arc_length = np.cos(lat1*PI_180) *2*np.pi*Re

    grid_lat_arc_length = np.sum(dy_[:,Ni//4]) 
    grid_lon_arc_length = np.sum(dx_[0,:])
    if(lat1>lat2):
        grid_lon_arc_length = np.sum(dx_[-1,:])
        
    area_error = 100*(np.sum(area_)-exact_area)/exact_area
    lat_arc_error = 100*(grid_lat_arc_length - exact_lat_arc_length)/exact_lat_arc_length
    lon_arc_error = 100*(grid_lon_arc_length -  exact_lon_arc_length)/exact_lon_arc_length
    print(exact_area)
    print(np.sum(area_))
    return lat_arc_error,lon_arc_error,area_error


def write_nc(x,y,dx,dy,area,angle_dx,axis_units='degrees',fnam=None,format='NETCDF3_64BIT',description=None,history=None,source=None,no_changing_meta=None):
    import netCDF4 as nc

    if fnam is None:
      fnam='supergrid.nc'
    fout=nc.Dataset(fnam,'w',format=format)

    ny=area.shape[0]; nx = area.shape[1]
    nyp=ny+1; nxp=nx+1
    print ('Writing netcdf file with ny,nx= ',ny,nx)

    nyp=fout.createDimension('nyp',nyp)
    nxp=fout.createDimension('nxp',nxp)
    ny=fout.createDimension('ny',ny)
    nx=fout.createDimension('nx',nx)
    string=fout.createDimension('string',255)    
    tile=fout.createVariable('tile','S1',('string'))
    yv=fout.createVariable('y','f8',('nyp','nxp'))
    xv=fout.createVariable('x','f8',('nyp','nxp'))    
    yv.units='degrees'
    xv.units='degrees'
    yv[:]=y
    xv[:]=x
#    tile[0:4]='tile1'
    tile[0]='t'
    tile[1]='i'
    tile[2]='l'
    tile[3]='e'
    tile[4]='1'
    dyv=fout.createVariable('dy','f8',('ny','nxp'))
    dyv.units='meters'
    dyv[:]=dy
    dxv=fout.createVariable('dx','f8',('nyp','nx'))
    dxv.units='meters'
    dxv[:]=dx
    areav=fout.createVariable('area','f8',('ny','nx'))
    areav.units='m2'
    areav[:]=area
    anglev=fout.createVariable('angle_dx','f8',('nyp','nxp'))
    anglev.units='degrees'
    anglev[:]=angle_dx
    #global attributes
    if(not no_changing_meta):
    	fout.history = history
    	fout.description = description
    	fout.source =  source

    fout.sync()
    fout.close()


def generate_latlon_grid(lni,lnj,llon0,llen_lon,llat0,llen_lat, ensure_nj_even=True):
    print('Generating regular lat-lon grid between latitudes ', llat0, llat0+llen_lat)
    llonSP = llon0 + np.arange(lni+1) * llen_lon/lni
    llatSP = llat0 + np.arange(lnj+1) * llen_lat/lnj
    if(llatSP.shape[0]%2 == 0 and ensure_nj_even):
        print("   The number of j's is not even. Fixing this by cutting one row at south.")
        llatSP = np.delete(llatSP,0,0)
    
    llamSP = np.tile(llonSP,(llatSP.shape[0],1)) 
    lphiSP = np.tile(llatSP.reshape((llatSP.shape[0],1)),(1,llonSP.shape[0]))
    
    print('   generated regular lat-lon grid between latitudes ', lphiSP[0,0],lphiSP[-1,0])
    print('   number of js=',lphiSP.shape[0])

    h_i_inv=llen_lon*PI_180*np.cos(lphiSP*PI_180)/lni
    h_j_inv=llen_lat*PI_180*np.ones(lphiSP.shape)/lnj
    delsin_j = np.roll(np.sin(lphiSP*PI_180),shift=-1,axis=0) - np.sin(lphiSP*PI_180)
    dx_h=h_i_inv[:,:-1]*_default_Re
    dy_h=h_j_inv[:-1,:]*_default_Re
    area=delsin_j[:-1,:-1]*_default_Re*_default_Re*llen_lon*PI_180/lni

    return llamSP,lphiSP,dx_h,dy_h,area

def usage():
    print('ocean_grid_generator.py -f <output_grid_filename> -r <inverse_degrees_resolution> [--rdp=<displacement_factor/0.2> --south_cutoff_ang=<degrees_south_to_start> --south_cutoff_row=<rows_south_to_cut> --reproduce_MIDAS_grids --reproduce_old8_grids --plot --write_subgrid_files]')
 

def main(argv):

    degree_resolution_inverse = 4 # (2 for half) or (4 for quarter) or (8 for 1/8) degree grid
    south_cap = True
    gridfilename = 'tripolar_res'+str(degree_resolution_inverse)+'.nc'
    doughnut=0.0
    r_dp=0.20
    rdp=1
    south_cutoff_row = 0
    south_cutoff_ang = -90.
    reproduce_MIDAS_grids = False
    reproduce_old8_grids = False
    write_subgrid_files = False
    plotem = False
    no_changing_meta = False

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hf:r:",["gridfilename=","inverse_resolution=","south_cutoff_ang=","south_cutoff_row=","rdp=","reproduce_MIDAS_grids","reproduce_old8_grids","plot","write_subgrid_files","no_changing_meta"])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)
        
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-f", "--gridfilename"):
            gridfilename = arg
        elif opt in ("-r", "--inverse_resolution"):
            degree_resolution_inverse = float(arg)
        elif opt in ("--south_cutoff_ang"):
            south_cutoff_ang = float(arg)
        elif opt in ("--south_cutoff_row"):
            south_cutoff_row = int(arg)
        elif opt in ("--rdp"):
             r_dp = int(arg)*r_dp
        elif opt in ("--reproduce_MIDAS_grids"):
             reproduce_MIDAS_grids = True
        elif opt in ("--reproduce_old8_grids"):
             reproduce_old8_grids = True
        elif opt in ("--plot"):
             plotem = True
        elif opt in ("--write_subgrid_files"):
             write_subgrid_files = True
        elif opt in ("--no_changing_meta"):
             no_changing_meta = True
        else:
            assert False, "unhandled option"


    #Information to write in file as metadata
    import socket
    host = str(socket.gethostname())
    scriptpath = sys.argv[0]
    scriptbasename = subprocess.check_output("basename "+ scriptpath,shell=True).decode('ascii').rstrip("\n")
    scriptdirname = subprocess.check_output("dirname "+ scriptpath,shell=True).decode('ascii').rstrip("\n")
    scriptgithash = subprocess.check_output("cd "+scriptdirname +";git rev-parse HEAD; exit 0",stderr=subprocess.STDOUT,shell=True).decode('ascii').rstrip("\n")
    scriptgitMod  = subprocess.check_output("cd "+scriptdirname +";git status --porcelain "+scriptbasename+" | awk '{print $1}' ; exit 0",stderr=subprocess.STDOUT,shell=True).decode('ascii').rstrip("\n")
    if("M" in str(scriptgitMod)):
        scriptgitMod = " , But was localy Modified!"

    hist = "This grid file was generated via command " + ' '.join(sys.argv)
    if(not no_changing_meta):
        hist = hist + " on "+ str(datetime.date.today()) + " on platform "+ host

    desc = "This is an orthogonal coordinate grid for the Earth with a nominal resoution of "+str(1/degree_resolution_inverse)+" degrees along the equator. "
    
    source =""
    if(not no_changing_meta):
        source =  source + scriptpath + " had git hash " + scriptgithash + scriptgitMod 
        source =  source + ". To obtain the grid generating code do: git clone  https://github.com/nikizadehgfdl/grid_generation.git ; cd grid_generation;  git checkout "+scriptgithash

    
    # Specify the grid properties
    # All
    # Specify the desired resolution
    refineS=2 # factor 2 is for supergrid
    refineR=degree_resolution_inverse   
    lenlon=360  # global longitude range
    lon0=-300.  # Starting longitude (longitude of the Northern bipoles)
    Ni=int(refineR*refineS* lenlon)
    #Apply finite element correction for cell areas when possible (currently for spherical sub-grids) 
    latlon_areafix=True
    #Ensure the number of j partitions are even for the sub-grids
    ensure_nj_even=True
    if(reproduce_old8_grids):
        latlon_areafix=False
 
    #MIDAS has nominal starting latitude for Mercator grid = -65 for 1/4 degree, -70 for 1/2 degree
    #MIDAS has nominal latitude range of Mercator grid     = 125 for 1/4 degree, 135 for 1/2 degree
    #Instead we use:
    phi_s_Merc, phi_n_Merc = -66.85954724706843, 64.0589597296948

    ###
    ###Mercator grid
    ###
    if(not reproduce_MIDAS_grids):
        lamMerc,phiMerc = generate_mercator_grid(Ni,phi_s_Merc,phi_n_Merc,lon0,lenlon, ensure_nj_even=ensure_nj_even)    
        dxMerc,dyMerc,areaMerc,angleMerc = generate_grid_metrics(lamMerc,phiMerc, latlon_areafix=latlon_areafix)        
    else: #use pymidas package   
        from pymidas.rectgrid_gen import supergrid

        if(refineR == 2):
            Nj_Merc    = 364*refineS  
            lat0_Merc  = -70.0  # This is a nominal starting latitude for Mercator grid
            lenlat_Merc= 135.0  # nominal latitude range of Mercator grid
        elif(refineR == 4):    
            Nj_Merc    = 700*refineS
            lat0_Merc  = -65.0  # This is a nominal starting latitude for Mercator grid
            lenlat_Merc= 125.0  # nominal latitude range of Mercator grid
        else:
            raise Exception('Unknown resolution for --reproduce_MIDAS_grids')
        
        #### Begin Mercator Grid
        print ('constructing a mercator supergrid with (ni,nj) = ',Ni,Nj_Merc)
        print ('nominal starting lat and starting longitude =',lat0_Merc, lon0)
        print ('and nominal width in latitude = ',lenlat_Merc)
        mercator=supergrid(Ni,Nj_Merc,'mercator','degrees',lat0_Merc,lenlat_Merc,lon0,360.,cyclic_x=True)

        #Add equatorial enhancement for 1/2 degree MIDAS grid
        if(refineR == 2):
            print ('Enhancing the equator resolution')
            import scipy.interpolate
            phi=np.ascontiguousarray( mercator.y[:,0] )
            dphi=phi[1:]-phi[0:-1]
            phi=mercator.y[:,0]
            dphi=phi[1:]-phi[0:-1]
            jind=np.where(phi>-30.)[0][0]
            jind=jind+np.mod(jind,2)
            phi=1.*phi[0:jind]
            dphi=dphi[0:jind]
            N=130
            phi_s = phi[-1]
            dphi_s = dphi[-1]
            phi_e = -5.
            dphi_e = 0.13
            nodes = [0,1,N-2,N-1]
            phi_nodes = [phi_s,phi_s+dphi_s,phi_e-dphi_e,phi_e]
            f2=scipy.interpolate.interp1d(nodes,phi_nodes,kind='cubic')
            jInd2=np.arange(N, dtype=float)
            phi2=f2(jInd2)
            print("Meridional range of pure Mercator=(", phi[0],",", phi[-2],") U (", -phi[-2],",", -phi[0],")." )
            print("Meridional range of cubic interpolation=(", phi2[0],"," , phi2[-2],") U (",-phi2[-2],",",-phi2[0],")." )
            phi=np.concatenate((phi[0:-1],phi2))
            N=40
            phi_s = phi[-1]
            phi2=np.linspace(phi_s,0,N)
            print("Meridional range of enhanced resolution=(", phi2[0],",", -phi2[0],").")
            print("Meridional value of enhanced resolution=", phi2[1]-phi2[0])
            PHI=np.concatenate((phi[0:-1],phi2))
            PHI=np.concatenate((PHI[0:-1],-PHI[::-1]))
            LAMBDA=np.linspace(lon0,lon0+360.,Ni+1)
            jind=np.where(PHI>-78.)[0][0]
            jind=jind+np.mod(jind,2)
            jind2=np.where(PHI>65.)[0][0]
            jind2=jind2+np.mod(jind2,2)
            PHI2=PHI[jind:jind2-1]
            x,y = np.meshgrid(LAMBDA,PHI2)
            mercator = supergrid(xdat=x,ydat=y,axis_units='degrees',cyclic_x=True)            

        print ("mercator.y.shape= ",mercator.y.shape)
        print ("mercator max/min latitude=", mercator.y.max(),mercator.y.min())
        print ("mercator start/end longitude=",mercator.x[0,0],mercator.x[0,-1])
        print ("mercator start/end latitudes=",mercator.y[0,0],mercator.y[-1,0])       
        mercator.grid_metrics()
        lamMerc,phiMerc = mercator.x,mercator.y
        dxMerc,dyMerc,areaMerc,angleMerc =mercator.dx,mercator.dy,mercator.area,mercator.angle_dx
            
    print("   CHECK_M: % errors in (lat arc, lon arc, area)", metrics_error(dxMerc,dyMerc,areaMerc,Ni,phiMerc[0,0],phiMerc[-1,0]))

    if(write_subgrid_files):
        write_nc(lamMerc,phiMerc,dxMerc,dyMerc,areaMerc,angleMerc,axis_units='degrees',fnam=gridfilename+"Merc.nc",description=desc,history=hist,source=source)

    #The phi resolution in the first and last row of Mercator grid along the symmetry meridian
    DeltaPhiMerc_so = phiMerc[ 1,Ni//4]-phiMerc[ 0,Ni//4]
    DeltaPhiMerc_no = phiMerc[-1,Ni//4]-phiMerc[-2,Ni//4]

    ###
    ###Northern bipolar cap
    ###
    lon_bp=lon0 # longitude of the displaced pole(s)

    if(not reproduce_MIDAS_grids):
        #Start lattitude from dy above the last Mercator grid
        lat0_bp = phiMerc[-1,Ni//4] + DeltaPhiMerc_no
        #Determine the number of bipolar cap grid point in the y direction such that the y resolution
        #along symmetry meridian is a constant and is equal to (continuous with) the last Mercator dy.
        #Note that int(0.5+x) is used to return the nearest integer to a float with deterministic behavior for middle points.
        #Note that int(0.5+x) is equivalent to math.floor(0.5+x)
        Nj_ncap = int(0.5+ (90.-lat0_bp)/DeltaPhiMerc_no) #Impose boundary condition for smooth dy
        
        if(reproduce_old8_grids):
            Nj_ncap=int(refineR* 120)
            lat0_bp=phi_n_Merc
        #Generate the bipolar grid
        lamBP,phiBP,dxBP_h,dyBP_h = generate_bipolar_cap_grid(Ni,Nj_ncap,lat0_bp,lon_bp,lenlon)
        #Metrics via MIDAS method
        dxBP,dyBP,areaBP,angleBP = generate_grid_metrics(lamBP,phiBP)
        print("   CHECK_M: % errors in (lat arc, lon arc, area)", metrics_error(dxBP,dyBP,areaBP,Ni,lat0_bp,90.))
        #Metrics via h's 
        dxBP_h=dxBP_h[:,:-1]
        dyBP_h=dyBP_h[:-1,:]
        arBP_h = dxBP_h[:-1,:]*dyBP_h[:,:-1]
        print("   CHECK_h0: % errors in (lat arc, lon arc, area)", metrics_error(dxBP_h,dyBP_h,arBP_h,Ni,lat0_bp,90.,Re=1))
        def quad_metrics(rf):
            #Metrics via h's, finite element quadrature
            lamBP_f,phiBP_f,dxBP_h_f,dyBP_h_f = generate_bipolar_cap_grid(Ni*rf,(Nj_ncap+1)*rf-1,lat0_bp,lon_bp,lenlon)
            print(lamBP_f.shape,phiBP_f.shape)
            print(dxBP_h_f.shape,dyBP_h_f.shape)
            dxBP_h_f=dxBP_h_f[:,:-1]        
            m,n=dxBP_h_f.shape
            dxBP_h_c=dxBP_h_f.reshape(m//rf,rf,n//rf,rf).sum(axis=(1,3))/rf

            dyBP_h_f=dyBP_h_f[:,:-1]
            m,n=dyBP_h_f.shape
            dyBP_h_c=dyBP_h_f.reshape(m//rf,rf,n//rf,rf).sum(axis=(1,3))/rf

            arBP_h_f = dxBP_h_f[:,:]*dyBP_h_f[:,:]
            arBP_h_c = arBP_h_f.reshape(m//rf,rf,n//rf,rf).sum(axis=(1,3))
            print("   CHECK_hq: % errors in (lat arc, lon arc, area)", metrics_error(dxBP_h_c,dyBP_h_c,arBP_h_c,Ni,lat0_bp,90.,Re=1))
        quad_metrics(10)
    else:
        if(refineR == 2):
            Nj_ncap = 119*refineS  
        elif(refineR == 4):    
            Nj_ncap = 240*refineS 
        else:
            raise Exception('Unknown resolution for --reproduce_MIDAS_grids')
            
        lat0_bp=mercator.y.max()
        dlat=90.0-lat0_bp

        tripolar_n=supergrid(Ni,Nj_ncap,'spherical','degrees',lat0_bp,dlat,lon0,360.,tripolar_n=True)
        tripolar_n.grid_metrics()
        print ("generated a tripolar supergrid of size (ny,nx)= ",tripolar_n.y.shape[0]-1,tripolar_n.y.shape[1]-1)
        print ("tripolar grid starting longitude = ",tripolar_n.x[0,0])
        print ("tripolar grid starting latitude = ",tripolar_n.y[0,0])

        lamBP,phiBP = tripolar_n.x,tripolar_n.y
        dxBP,dyBP,areaBP,angleBP =tripolar_n.dx,tripolar_n.dy,tripolar_n.area,tripolar_n.angle_dx


    if(write_subgrid_files):
        write_nc(lamBP,phiBP,dxBP,dyBP,areaBP,angleBP,axis_units='degrees',fnam=gridfilename+"BP.nc",description=desc,history=hist,source=source)
    ###
    ###Southern Ocean grid
    ###
    lat0_SO=-78. #Starting lat

    if(not reproduce_MIDAS_grids):
        #Make the last grid point is a (Mercator) step below the first Mercator lattitude.
        lenlat_SO = phiMerc[0,Ni//4] - DeltaPhiMerc_so - lat0_SO #Start from a lattitude to smooth out dy.
        #Determine the number of grid point in the y direction such that the y resolution is equal to (continuous with)
        #the first Mercator dy.     
        Nj_SO = int(0.5 + lenlat_SO/DeltaPhiMerc_so) #Make the resolution continious with the Mercator at joint

        if(reproduce_old8_grids):
            lenlat_SO = phi_s_Merc-lat0_SO
            Nj_SO  =int(refineR*  55)

        lamSO,phiSO,dxhSO,dyhSO,arhSO = generate_latlon_grid(Ni,Nj_SO,lon0,lenlon,lat0_SO,lenlat_SO, ensure_nj_even=ensure_nj_even)
        dxSO,dySO,areaSO,angleSO = generate_grid_metrics(lamSO,phiSO, latlon_areafix=latlon_areafix)
        #Metrics errors via MIDAS
        print("   CHECK_M: % errors in (lat arc, lon arc, area)", metrics_error(dxSO,dySO,areaSO,Ni,phiSO[0,0],phiSO[-1,0]))
        #Metrics errors via h's 
        print("   CHECK_h: % errors in (lat arc, lon arc, area)", metrics_error(dxhSO,dyhSO,arhSO,Ni,phiSO[0,0],phiSO[-1,0]))
        
    else:
        if(refineR == 2):
            Nj_SO = 54*refineS  
        elif(refineR == 4):    
            Nj_SO = 110*refineS 
        else:
            raise Exception('Unknown resolution for --reproduce_MIDAS_grids. Use either -r 2 or -r 4 ')
            
        print ('constructing a spherical supergrid with (ny,nx) = ',Ni,Nj_SO)
        print ('nominal starting lat and starting longitude =',lat0_SO, lon0)
        print ('and nominal width in latitude = ',mercator.y.min()-lat0_SO)
        spherical=supergrid(Ni,Nj_SO,'spherical','degrees',lat0_SO,mercator.y.min()-lat0_SO,lon0,360.,cyclic_x=True)
        spherical.grid_metrics()
        print ("southern ocean spherical max/min latitude=", spherical.y.max(),spherical.y.min())
        print ("spherical nj,ni=", spherical.y.shape[0]-1,spherical.y.shape[1]-1)
        print ("spherical starting longitude=",spherical.x[0,0])
        print ("spherical ending longitude=",spherical.x[0,-1])

        lamSO,phiSO = spherical.x,spherical.y
        dxSO,dySO,areaSO,angleSO =spherical.dx,spherical.dy,spherical.area,spherical.angle_dx

    if(write_subgrid_files):
        write_nc(lamSO,phiSO,dxSO,dySO,areaSO,angleSO,axis_units='degrees',fnam=gridfilename+"SO.nc",description=desc,history=hist,source=source)
    ###
    ###Southern cap
    ###
    #Nj_scap = refineR *  40   #MIDAS has refineS*(  80 for 1/4 degree, ??? for 1/2 degree
    #Here we determine Nj_scap by imposing the condition of continuity of dy across the stitch.
    
    lon_dp=100.0   # longitude of the displaced pole 
    deltaPhiSO = phiSO[1,Ni//4]-phiSO[0,Ni//4]
    lat0_SC=phiSO[0,Ni//4]-deltaPhiSO
    fullArc = lat0_SC+90.
    halfArc = fullArc/2

    if(r_dp == 0.0):
        if(not reproduce_MIDAS_grids):
            Nj_scap = int(fullArc/deltaPhiSO)
            if(reproduce_old8_grids):
                Nj_scap=int(refineR*  40)
                lat0_SC=lat0_SO

            lamSC,phiSC,dxhSC,dyhSC,arhSC = generate_latlon_grid(Ni,Nj_scap,lon0,lenlon,-90.,90+lat0_SO, ensure_nj_even=ensure_nj_even)
            dxSC,dySC,areaSC,angleSC = generate_grid_metrics(lamSC,phiSC, latlon_areafix=latlon_areafix)
            #Metrics errors via h's 
            print("   CHECK_h: % errors in (lat arc, lon arc, area)", metrics_error(dxhSC,dyhSC,arhSC,Ni,phiSC[-1,0],phiSC[0,0]))
        else:
            Nj_scap = int(refineR *  40)   #MIDAS has refineS*(  80 for 1/4 degree, ??? for 1/2 degree
            spherical_cap=supergrid(Ni,Nj_scap,'spherical','degrees',-90.,lat0_SO+90,lon0,360.,cyclic_x=True)
            spherical_cap.grid_metrics()
            print ("spherical cap max/min latitude=", spherical_cap.y.max(),spherical_cap.y.min())
            print ("spherical cap nj,ni=", spherical_cap.y.shape[0]-1,spherical_cap.y.shape[1]-1)
            print ("spherical cap starting longitude=",spherical_cap.x[0,0])
            print ("spherical cap ending longitude=",spherical_cap.x[0,-1])

            lamSC,phiSC = spherical_cap.x,spherical_cap.y
            dxSC,dySC,areaSC,angleSC =spherical_cap.dx,spherical_cap.dy,spherical_cap.area,spherical_cap.angle_dx

    else:
        doughnut=0.12
        nparts=8
        Nj_scap = int((nparts/(nparts-1))*halfArc/deltaPhiSO)
        if(not reproduce_MIDAS_grids):
            if(reproduce_old8_grids):
                Nj_scap=int(refineR*  40)
                lat0_SC=lat0_SO

            lamSC,phiSC = generate_displaced_pole_grid(Ni,Nj_scap,lon0,lenlon,lon_dp,r_dp,lat0_SC,doughnut,nparts, ensure_nj_even=ensure_nj_even)
            dxSC,dySC,areaSC,angleSC = generate_grid_metrics(lamSC,phiSC)
            
        else:    
            Nj_scap = int(refineR *  40)   #MIDAS has refineS*(  80 for 1/4 degree, ??? for 1/2 degree
            ny_scap = Nj_scap
            r0_pole = 0.20
            lon0_pole=100.0

            print( spherical.dy.shape)
            lenlat=90.0+spherical.y.min()
            dy0=spherical.dy[0,0]*r0_pole
            x=spherical.x[0,:]
            y=np.linspace(-90.,0.5*(lat0_SO-90.0),ny_scap/8)
            y=np.concatenate((y,np.linspace(y.max(),lat0_SO,7*ny_scap/8+1)))
            X,Y=np.meshgrid(x,y)
            antarctic_cap=supergrid(xdat=X,ydat=Y,axis_units='degrees',displace_pole=True,r0_pole=r0_pole,lon0_pole=lon0_pole,doughnut=doughnut)                  
            antarctic_cap.grid_metrics()
            print( "generated a southern cap of size (ny,nx)= ",antarctic_cap.y.shape[0]-1,antarctic_cap.y.shape[1]-1)

            lamSC,phiSC = antarctic_cap.x,antarctic_cap.y
            dxSC,dySC,areaSC,angleSC =antarctic_cap.dx,antarctic_cap.dy,antarctic_cap.area,antarctic_cap.angle_dx

    print("   CHECK_M: % errors in (lat arc, lon arc, area)", metrics_error(dxSC,dySC,areaSC,Ni,phiSC[-1,0],phiSC[0,0]))
    if(write_subgrid_files):
        write_nc(lamSC,phiSC,dxSC,dySC,areaSC,angleSC,axis_units='degrees',fnam=gridfilename+"SC.nc",description=desc,history=hist,source=source)
    #Concatenate to generate the whole grid
    #Start from displaced southern cap and join the southern ocean grid
    print("Stitching the grids together...")
    x1=np.concatenate((lamSC,lamSO[1:,:]),axis=0)
    y1=np.concatenate((phiSC,phiSO[1:,:]),axis=0)
    dx1=np.concatenate((dxSC,dxSO[1:,:]),axis=0)
    dy1=np.concatenate((dySC,dySO),axis=0)
    area1=np.concatenate((areaSC,areaSO),axis=0)
    angle1=np.concatenate((angleSC[:-1,:],angleSO[:-1,:]),axis=0)
    #Join the Mercator grid
    x2=np.concatenate((x1,lamMerc[1:,:]),axis=0)
    y2=np.concatenate((y1,phiMerc[1:,:]),axis=0)
    dx2=np.concatenate((dx1,dxMerc[1:,:]),axis=0)
    dy2=np.concatenate((dy1,dyMerc),axis=0)
    area2=np.concatenate((area1,areaMerc),axis=0)
    angle2=np.concatenate((angle1,angleMerc[:-1,:]),axis=0)
    #Join the norhern bipolar cap grid
    x3=np.concatenate((x2,lamBP[1:,:]),axis=0)
    y3=np.concatenate((y2,phiBP[1:,:]),axis=0)
    dx3=np.concatenate((dx2,dxBP[1:,:]),axis=0)
    dy3=np.concatenate((dy2,dyBP),axis=0)
    area3=np.concatenate((area2,areaBP),axis=0)
    angle3=np.concatenate((angle2,angleBP),axis=0)

    ##Drop the first 80 points from the southern cap! Make this an option!
    if(south_cutoff_row > 0):
        jcut = south_cutoff_row - 1
        print("cutting grid rows 0 to ", jcut)
        x3 = x3[jcut:,:]
        y3 = y3[jcut:,:]
        dx3 = dx3[jcut:,:]
        dy3 = dy3[jcut:,:]
        area3 = area3[jcut:,:]
        angle3 = angle3[jcut:,:]
        
    if(south_cutoff_ang > -90):
        jcut = 1 + np.nonzero(y3[:,0] < south_cutoff_ang)[0][-1]
        print("cutting grid below ",south_cutoff_ang, jcut)
        x3 = x3[jcut:,:]
        y3 = y3[jcut:,:]
        dx3 = dx3[jcut:,:]
        dy3 = dy3[jcut:,:]
        area3 = area3[jcut:,:]
        angle3 = angle3[jcut:,:]
        
        
    #write the whole grid file    
    desc = desc + "It consists of a Mercator grid spanning "+ str(phiMerc[0,0]) + " to " + str(phiMerc[-1,0]) + " degrees, flanked by a bipolar northern cap and a regular lat-lon grid spanning " + str(phiMerc[0,0]) + " to " + str(lat0_SO)+" degrees. "

    desc = desc + "It is capped by a "
    if(r_dp != 0.0): 
        desc = desc + "displaced pole "
    else:
        desc = desc + "regular "
    desc = desc + "grid south of "+ str(lat0_SO) +" degrees."
    
    if(south_cutoff_ang > -90):
        desc = desc + " It is cut south of " + str(south_cutoff_ang) + " degrees."

    if(south_cutoff_row > 0):
        desc = desc + " The first "+ str(south_cutoff_row) +" rows at south are deleted."
    

#    print(hist)
#    print(desc)
#    print(source)

    write_nc(x3,y3,dx3,dy3,area3,angle3,axis_units='degrees',fnam=gridfilename,description=desc,history=hist,source=source,no_changing_meta=no_changing_meta)
    print("Wrote the whole grid to file ",gridfilename)
    
    #Visualization
    if(plotem):
        plot_mesh_in_xyz(x2,y2, stride=30,upperlat=-40, title="Grid south of -40 degrees")
        plot_mesh_in_xyz(x3,y3, stride=30,lowerlat=40, title="Grid north of 40 degrees") 

if __name__ == "__main__":
   main(sys.argv[1:])

