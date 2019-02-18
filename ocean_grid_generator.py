#!/usr/bin/env python

import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import numpy as np
import sys, getopt
import datetime, os, subprocess

#Constants
PI_180 = np.pi/180.
#_default_Re = 6.378e6
_default_Re = 6371.e3 #MIDAS
#_default_Re = 6.371e6

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
    return lams,phis

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
            print( 'y*=',y_star, 'nj=', y_star[1]-y_star[0]+1 )
    if((y_star[1]-y_star[0]+1)%2 != 0):
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
    if(y.shape[0]%2 != 0 and ensure_nj_even):
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
#plot_mesh_in_latlon(lams,phis,stride=16)

def plot_mesh_in_xyz(lam, phi, stride=1, phi_color='k', lam_color='r', lowerlat=None, upperlat=None, newfig=True, title=None):
    if lowerlat is not None:
        lam,phi = cut_below(lam,phi,lowerlat=lowerlat)
        
    if upperlat is not None:
        lam,phi = cut_above(lam,phi,upperlat=upperlat)
        
    x = np.cos(phi*PI_180) * np.cos(lam*PI_180)
    y = np.cos(phi*PI_180) * np.sin(lam*PI_180)
    z = np.sin(phi*PI_180)
        
    plot_mesh_in_latlon(x, y, stride=stride, phi_color=phi_color, lam_color=lam_color, newfig=newfig, title=title)
#plt.figure(figsize=(6,6))
#plot_mesh_in_xyz(lams, phis, stride=20)

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


def write_nc(x,y,dx,dy,area,angle_dx,axis_units='degrees',fnam=None,format='NETCDF3_CLASSIC',description=None,history=None,source=None):
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
    fout.history = history
    fout.description = description
    fout.source =  source

    fout.sync()
    fout.close()


def generate_latlon_grid(lni,lnj,llon0,llen_lon,llat0,llen_lat, ensure_nj_even=True):
    print('Generating regular lat-lon grid between latitudes ', llat0, llat0+llen_lat)
    llonSP = llon0 + np.arange(lni+1) * llen_lon/lni
    llatSP = llat0 + np.arange(lnj+1) * llen_lat/lnj
    if(llatSP.shape[0]%2 != 0 and ensure_nj_even):
        print("   The number of j's is not even. Fixing this by cutting one row at south.")
        llatSP = np.delete(llatSP,0,0)
    
    llamSP = np.tile(llonSP,(llatSP.shape[0],1)) 
    lphiSP = np.tile(llatSP.reshape((llatSP.shape[0],1)),(1,llonSP.shape[0]))
    
    print('   generated regular lat-lon grid between latitudes ', lphiSP[0,0],lphiSP[-1,0])
    print('   number of js=',lphiSP.shape[0])
    return llamSP,lphiSP

def usage():
    print('ocean_grid_generator.py -f <output_grid_filename> -r <inverse_degrees_resolution> --rdp=<displacement_factor/0.2> --south_cutoff_ang=<degrees_south_to_start> --south_cutoff_row=<rows_south_to_cut>')
 

def main(argv):

    degree_resolution_inverse = 4 # (2 for half) or (4 for quarter) or (8 for 1/8) degree grid
    south_cap = True
    gridfilename = 'tripolar_res'+str(degree_resolution_inverse)+'.nc'
    doughnut=0.0
    r_dp=0.20
    rdp=1
    south_cutoff_row = 0
    south_cutoff_ang = -90.
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hf:r:",["gridfilename=","inverse_resolution=","south_cutoff_ang=","south_cutoff_row=","rdp="])
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
            degree_resolution_inverse = int(arg)
        elif opt in ("--south_cutoff_ang"):
            south_cutoff_ang = float(arg)
        elif opt in ("--south_cutoff_row"):
            south_cutoff_row = int(arg)
        elif opt in ("--rdp"):
             r_dp = int(arg)*r_dp
        else:
            assert False, "unhandled option"

            
    scriptpath = sys.argv[0]
#    print("scriptpath=",scriptpath)   
    scriptbasename = subprocess.check_output("basename "+ scriptpath,shell=True).decode('ascii').rstrip("\n")
    scriptdirname = subprocess.check_output("dirname "+ scriptpath,shell=True).decode('ascii').rstrip("\n")
#    print("scriptname=",scriptname)
#    print(subprocess.check_output("/usr/local/x64/git-2.4.6/bin/git rev-parse HEAD; /usr/local/x64/git-2.4.6/bin/git status --porcelain `basename "+sys.argv[0]+"` ; exit 0",stderr=subprocess.STDOUT,shell=True))
#    scriptgithash = subprocess.check_output("/usr/local/x64/git-2.4.6/bin/git rev-parse HEAD; /usr/local/x64/git-2.4.6/bin/git status --porcelain `basename "+sys.argv[0]+"` ; exit 0",stderr=subprocess.STDOUT,shell=True)
    scriptgithash = subprocess.check_output("cd "+scriptdirname +";/usr/local/x64/git-2.4.6/bin/git rev-parse HEAD; exit 0",stderr=subprocess.STDOUT,shell=True).decode('ascii').rstrip("\n")
    scriptgitMod  = subprocess.check_output("/usr/local/x64/git-2.4.6/bin/git status --porcelain "+scriptbasename+" | awk '{print $1}' ; exit 0",stderr=subprocess.STDOUT,shell=True).decode('ascii').rstrip("\n")
    if("M" in str(scriptgitMod)):
        scriptgitMod = " , But was localy Modified!"
    
    # Specify the grid properties
    # All
    # Specify the desired resolution
    refineS=2 # factor 2 is for supergrid
    refineR=degree_resolution_inverse   
    lenlon=360  # global longitude range
    lon0=-300.  # Starting longitude (longitude of the Northern bipoles)
    Ni     =refineR*refineS* lenlon

    #Mercator grid
    #MIDAS has nominal starting latitude for Mercator grid = -65 for 1/4 degree, -70 for 1/2 degree
    #MIDAS has nominal latitude range of Mercator grid     = 125 for 1/4 degree, 135 for 1/2 degree
    #Instead we use:
    phi_s_Merc, phi_n_Merc = -66.85954724706843, 64.0589597296948

#    if(degree_resolution_inverse == 2):
#        phi_s_Merc, phi_n_Merc = -68.05725377, 64.76080366
#    OM4p5 grid is not Mercator around equator and cannot be reproduced by this tool
     
    lamMerc,phiMerc = generate_mercator_grid(Ni,phi_s_Merc,phi_n_Merc,lon0,lenlon,ensure_nj_even=True)    
    dxMerc,dyMerc,areaMerc,angleMerc = generate_grid_metrics(lamMerc,phiMerc,axis_units='degrees',latlon_areafix=True)
    #The phi resolution in the first and last row of Mercator grid along the symmetry meridian
    DeltaPhiMerc_so = phiMerc[ 1,Ni//4]-phiMerc[ 0,Ni//4]
    DeltaPhiMerc_no = phiMerc[-1,Ni//4]-phiMerc[-2,Ni//4]

    #Northern bipolar cap
    lon_bp=lon0 # longitude of the displaced pole(s)
    #Start lattitude from dy above the last Mercator grid
    lat0_bp = phiMerc[-1,Ni//4] + DeltaPhiMerc_no
    #Determine the number of bipolar cap grid point in the y direction such that the y resolution
    #along symmetry meridian is a constant and is equal to (continuous with) the last Mercator dy.
    #Note that int(0.5+x) is used to return the nearest integer to a float with deterministic behavior for middle points.
    #Note that int(0.5+x) is equivalent to math.floor(0.5+x)
    Nj_ncap = int(0.5+ (90.-lat0_bp)/DeltaPhiMerc_no) #Impose boundary condition for smooth dy

    lamBP,phiBP = generate_bipolar_cap_grid(Ni,Nj_ncap,lat0_bp,lon_bp,lenlon)
    dxBP,dyBP,areaBP,angleBP = generate_grid_metrics(lamBP,phiBP,axis_units='degrees')

    #Southern Ocean grid
    lat0_SO=-78. #Starting lat
    #Make the last grid point is a (Mercator) step below the first Mercator lattitude.
    lenlat_SO = phiMerc[0,Ni//4] - DeltaPhiMerc_so - lat0_SO #Start from a lattitude to smooth out dy.
    #Determine the number of grid point in the y direction such that the y resolution is equal to (continuous with)
    #the first Mercator dy.     
    Nj_SO = int(0.5 + lenlat_SO/DeltaPhiMerc_so) #Make the resolution continious with the Mercator at joint
    lamSO,phiSO = generate_latlon_grid(Ni,Nj_SO,lon0,lenlon,lat0_SO,lenlat_SO)
    dxSO,dySO,areaSO,angleSO = generate_grid_metrics(lamSO,phiSO,axis_units='degrees',latlon_areafix=True)

    #Southern cap
    #Nj_scap = refineR *  40   #MIDAS has refineS*(  80 for 1/4 degree, ??? for 1/2 degree
    #Here we determine Nj_scap by imposing the condition of continuity of dy across the stitch.
    
    lon_dp=100.0   # longitude of the displaced pole 
    deltaPhiSO = phiSO[1,Ni//4]-phiSO[0,Ni//4]
    lat0_SC=phiSO[0,Ni//4]-deltaPhiSO
    fullArc = lat0_SC+90.
    halfArc = fullArc/2

    if(r_dp == 0.0):
        Nj_scap = int(fullArc/deltaPhiSO)
        lamSC,phiSC = generate_latlon_grid(Ni,Nj_scap,lon0,lenlon,-90.,90+lat0_SO)
    else:
        nparts=8
        Nj_scap = int((nparts/(nparts-1))*halfArc/deltaPhiSO)
        lamSC,phiSC = generate_displaced_pole_grid(Ni,Nj_scap,lon0,lenlon,lon_dp,r_dp,lat0_SC,doughnut,nparts)

    dxSC,dySC,areaSC,angleSC = generate_grid_metrics(lamSC,phiSC,axis_units='degrees')

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
        
        
        
    #write the grid file    
    hist = "This grid file was generated on "+ str(datetime.date.today()) + " via command " + ' '.join(sys.argv)
    desc = "This is an orthogonal coordinate grid for the Earth with a nominal resoution of "+str(1/degree_resolution_inverse)+" degrees along the equator. "

    desc = desc + "It consists of a Mercator grid spanning "+ str(phi_s_Merc) + " to " + str(phi_n_Merc) + " degrees, flanked by a bipolar northern cap and a regular lat-lon grid spanning " + str(phi_s_Merc) + " to " + str(lat0_SO)+" degrees. "

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
    
    source =  scriptpath + " had git hash " + scriptgithash + scriptgitMod 
    source =  source + ". To obtain the grid generating code do: git clone  https://github.com/nikizadehgfdl/grid_generation.git ; cd grid_generation;  git checkout "+scriptgithash

#    print(hist)
#    print(desc)
#    print(source)

    write_nc(x3,y3,dx3,dy3,area3,angle3,axis_units='degrees',fnam=gridfilename,description=desc,history=hist,source=source)
    print("Wrote the whole grid to file ",gridfilename)
    

    #Visualization
#    plot_mesh_in_xyz(x2,y2, stride=30,upperlat=-40, title="Grid south of -40 degrees")
#    plot_mesh_in_xyz(x3,y3, stride=30,lowerlat=40, title="Grid north of 40 degrees")
#    plt.show()

    


if __name__ == "__main__":
   main(sys.argv[1:])
