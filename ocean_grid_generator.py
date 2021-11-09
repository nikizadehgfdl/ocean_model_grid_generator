#!/usr/bin/env python

from __future__ import print_function

import numpypi.numpypi_series as np

# import numpy as np
import sys, getopt
import datetime, os, subprocess

# Constants
PI_180 = np.pi / 180.0
# _default_Re = 6.378e6
_default_Re = 6371.0e3  # MIDAS
HUGE = 1.0e30


def chksum(x, lbl):
    import hashlib

    if type(x) in (float, int, np.float64):
        y = np.array(x)
    else:
        y = np.zeros(x.shape)
        y[:] = x
    ymin, ymax, ymean = y.min(), y.max(), y.mean()
    ysd = np.sqrt(((y - ymean) ** 2).mean())
    print(
        hashlib.sha256(y).hexdigest(),
        "%10s" % lbl,
        "min = %.15f" % ymin,
        "max = %.15f" % ymax,
        "mean = %.15f" % ymean,
        "sd = %.15f" % ysd,
    )


def bipolar_projection(lamg, phig, lon_bp, rp, metrics_only=False):
    """Makes a stereographic bipolar projection of the input coordinate mesh (lamg,phig)
    Returns the projected coordinate mesh and their metric coefficients (h^-1).
    The input mesh must be a regular spherical grid capping the pole with:
        latitudes between 2*arctan(rp) and 90  degrees
        longitude between lon_bp       and lonp+360
    """
    ### symmetry meridian resolution fix
    phig = 90 - 2 * np.arctan(np.tan(0.5 * (90 - phig) * PI_180) / rp) / PI_180
    tmp = mdist(lamg, lon_bp) * PI_180
    sinla = np.sin(tmp)  # This makes phis symmetric
    sphig = np.sin(phig * PI_180)
    alpha2 = (np.cos(tmp)) ** 2  # This makes dy symmetric
    beta2_inv = (np.tan(phig * PI_180)) ** 2
    rden = 1.0 / (1.0 + alpha2 * beta2_inv)

    if not metrics_only:
        B = sinla * np.sqrt(rden)  # Actually two equations  +- |B|
        # Deal with beta=0
        B = np.where(np.abs(beta2_inv) > HUGE, 0.0, B)
        lamc = np.arcsin(B) / PI_180
        ##But this equation accepts 4 solutions for a given B, {l, 180-l, l+180, 360-l }
        ##We have to pickup the "correct" root.
        ##One way is simply to demand lamc to be continuous with lam on the equator phi=0
        ##I am sure there is a more mathematically concrete way to do this.
        lamc = np.where((lamg - lon_bp > 90) & (lamg - lon_bp <= 180), 180 - lamc, lamc)
        lamc = np.where(
            (lamg - lon_bp > 180) & (lamg - lon_bp <= 270), 180 + lamc, lamc
        )
        lamc = np.where((lamg - lon_bp > 270), 360 - lamc, lamc)
        # Along symmetry meridian choose lamc
        lamc = np.where(
            (lamg - lon_bp == 90), 90, lamc
        )  # Along symmetry meridian choose lamc=90-lon_bp
        lamc = np.where(
            (lamg - lon_bp == 270), 270, lamc
        )  # Along symmetry meridian choose lamc=270-lon_bp
        lams = lamc + lon_bp

    ##Project back onto the larger (true) sphere so that the projected equator shrinks to latitude \phi_P=lat0_tp
    ##then we have tan(\phi_s'/2)=tan(\phi_p'/2)tan(\phi_c'/2)
    A = sinla * sphig
    chic = np.arccos(A)
    phis = 90 - 2 * np.arctan(rp * np.tan(chic / 2)) / PI_180
    ##Calculate the Metrics
    rden2 = 1.0 / (1 + (rp * np.tan(chic / 2)) ** 2)
    M_inv = rp * (1 + (np.tan(chic / 2)) ** 2) * rden2
    M = 1 / M_inv
    chig = (90 - phig) * PI_180
    rden2 = 1.0 / (1 + (rp * np.tan(chig / 2)) ** 2)
    N = rp * (1 + (np.tan(chig / 2)) ** 2) * rden2
    N_inv = 1 / N
    cos2phis = (np.cos(phis * PI_180)) ** 2

    h_j_inv = (
        cos2phis * alpha2 * (1 - alpha2) * beta2_inv * (1 + beta2_inv) * (rden ** 2)
        + M_inv * M_inv * (1 - alpha2) * rden
    )
    # Deal with beta=0. Prove that cos2phis/alpha2 ---> 0 when alpha, beta  ---> 0
    h_j_inv = np.where(np.abs(beta2_inv) > HUGE, M_inv * M_inv, h_j_inv)
    h_j_inv = np.sqrt(h_j_inv) * N_inv

    h_i_inv = (
        cos2phis * (1 + beta2_inv) * (rden ** 2)
        + M_inv * M_inv * alpha2 * beta2_inv * rden
    )
    # Deal with beta=0
    h_i_inv = np.where(np.abs(beta2_inv) > HUGE, M_inv * M_inv, h_i_inv)
    h_i_inv = np.sqrt(h_i_inv)

    if not metrics_only:
        return lams, phis, h_i_inv, h_j_inv
    else:
        return h_i_inv, h_j_inv


def generate_bipolar_cap_mesh(Ni, Nj_ncap, lat0_bp, lon_bp, ensure_nj_even=True):
    # Define a (lon,lat) coordinate mesh on the Northern hemisphere of the globe sphere
    # such that the resolution of latg matches the desired resolution of the final grid along the symmetry meridian
    print("Generating bipolar grid bounded at latitude ", lat0_bp)
    if Nj_ncap % 2 != 0 and ensure_nj_even:
        print("   Supergrid has an odd number of area cells!")
        if ensure_nj_even:
            print("   The number of j's is not even. Fixing this by cutting one row.")
            Nj_ncap = Nj_ncap - 1

    lon_g = lon_bp + np.arange(Ni + 1) * 360.0 / float(Ni)
    lamg = np.tile(lon_g, (Nj_ncap + 1, 1))
    latg0_cap = lat0_bp + np.arange(Nj_ncap + 1) * (90 - lat0_bp) / float(Nj_ncap)
    phig = np.tile(latg0_cap.reshape((Nj_ncap + 1, 1)), (1, Ni + 1))
    rp = np.tan(0.5 * (90 - lat0_bp) * PI_180)
    lams, phis, h_i_inv, h_j_inv = bipolar_projection(lamg, phig, lon_bp, rp)
    h_i_inv = h_i_inv[:, :-1] * 2 * np.pi / float(Ni)
    h_j_inv = h_j_inv[:-1, :] * PI_180 * (90 - lat0_bp) / float(Nj_ncap)
    print("   number of js=", phis.shape[0])
    return lams, phis, h_i_inv, h_j_inv


def bipolar_cap_ij_array(i, j, Ni, Nj_ncap, lat0_bp, lon_bp, rp):
    long = lon_bp + i * 360.0 / float(Ni)
    latg = lat0_bp + j * (90 - lat0_bp) / float(Nj_ncap)
    lamg = np.tile(long, (latg.shape[0], 1))
    phig = np.tile(latg.reshape((latg.shape[0], 1)), (1, long.shape[0]))
    h_i_inv, h_j_inv = bipolar_projection(lamg, phig, lon_bp, rp, metrics_only=True)
    h_i_inv = h_i_inv * 2 * np.pi / float(Ni)
    h_j_inv = h_j_inv * (90 - lat0_bp) * PI_180 / float(Nj_ncap)
    return h_i_inv, h_j_inv


def bipolar_cap_metrics_quad_fast(order, nx, ny, lat0_bp, lon_bp, rp, Re=_default_Re):
    print("   Calculating bipolar cap metrics via quadrature ...")
    a, b = quad_positions(order)
    daq = np.zeros([ny + 1, nx + 1])
    dxq = np.zeros([ny + 1, nx + 1])
    dyq = np.zeros([ny + 1, nx + 1])

    j1d = np.empty([0])
    for j in range(0, ny + 1):
        j_s = b * j + a * (j + 1)
        if j_s[-1] == ny:
            j_s[-1] = ny - 0.001  # avoid phi=90 as this will cause errore.
        # Niki:Find a way to avoid this properly.
        # This could be a sign that there is still something
        # wrong with the h_j_inv calculations at phi=90 (beta=0).
        j1d = np.append(j1d, j_s)

    i1d = np.empty([0])
    for i in range(0, nx + 1):
        i_s = b * i + a * (i + 1)
        i1d = np.append(i1d, i_s)

    # dx,dy = bipolar_cap_ij_array(i1d,j1d,nx,ny,lat0_bp,lon_bp,rp)
    # Or to make it faster:
    nj, ni = j1d.shape[0], i1d.shape[0]  # Shape of results
    dj = min(nj, max(32 * 1024 // ni, 1))  # Stride to use that fits in memory
    lams, phis, dx, dy = (
        np.zeros((nj, ni)),
        np.zeros((nj, ni)),
        np.zeros((nj, ni)),
        np.zeros((nj, ni)),
    )
    for j in range(0, nj, dj):
        je = min(nj, j + dj)
        dx[j:je], dy[j:je] = bipolar_cap_ij_array(
            i1d, j1d[j:je], nx, ny, lat0_bp, lon_bp, rp
        )

    # reshape to send for quad averaging
    dx_r = dx.reshape(ny + 1, order, nx + 1, order)
    dy_r = dy.reshape(ny + 1, order, nx + 1, order)
    # area element
    dxdy_r = dx_r * dy_r

    for j in range(0, ny + 1):
        for i in range(0, nx + 1):
            daq[j, i] = quad_average_2d(dxdy_r[j, :, i, :])
            dxq[j, i] = quad_average(dx_r[j, 0, i, :])
            dyq[j, i] = quad_average(dy_r[j, :, i, 0])
    daq = daq[:-1, :-1] * Re * Re
    dxq = dxq[:, :-1] * Re
    dyq = dyq[:-1, :] * Re
    return dxq, dyq, daq


def quad_positions(n=3):
    """Returns weights wa and wb so that the element [xa,xb] is sampled at positions
    x=wa(xa+xb*xb)."""
    if n == 2:
        return np.array([0.0, 1.0]), np.array([1.0, 0.0])
    if n == 3:
        return np.array([0.0, 0.5, 1.0]), np.array([1.0, 0.5, 0.0])
    if n == 4:
        r5 = 0.5 / np.sqrt(5.0)
        return np.array([0.0, 0.5 - r5, 0.5 + r5, 1.0]), np.array(
            [1.0, 0.5 + r5, 0.5 - r5, 0.0]
        )
    if n == 5:
        r37 = 0.5 * np.sqrt(3.0 / 7.0)
        return np.array([0.0, 0.5 - r37, 0.5, 0.5 + r37, 1.0]), np.array(
            [1.0, 0.5 + r37, 0.5, 0.5 - r37, 0.0]
        )
    raise Exception("Uncoded order")


def quad_average(y):
    """Returns the average value found by quadrature at order n.
    y is a list of values in order from x=-1 to x=1."""
    if len(y) == 2:  # 1, 1
        d = 1.0 / 2.0
        return d * (y[0] + y[1])
    if len(y) == 3:  # 1/3, 4/3, 1/3
        d = 1.0 / 6.0
        return d * (4.0 * y[1] + (y[0] + y[2]))
    if len(y) == 4:  # 1/6, 5/6, 5/6, 1/6
        d = 1.0 / 12.0
        return d * (5.0 * (y[1] + y[2]) + (y[0] + y[3]))
    if len(y) == 5:  # 9/10, 49/90, 64/90, 49/90, 9/90
        d = 1.0 / 180.0
        return d * (64.0 * y[2] + (49.0 * (y[1] + y[3])) + 9.0 * (y[0] + y[4]))
    raise Exception("Uncoded order")


def quad_average_2d(y):
    """Returns the average value found by quadrature at order n.
    y is a list of values in order from x1=-1 to x1=1 and x2=-1 to x2=1."""
    if y.shape[0] != y.shape[1]:
        raise Exception("Input array is not squared!")

    if y.shape[0] == 2:  # 1, 1
        d = 1.0 / 2.0
        return d * d * (y[0, 0] + y[0, 1] + y[1, 0] + y[1, 1])
    if y.shape[0] == 3:  # 1/3, 4/3, 1/3
        d = 1.0 / 6.0
        return (
            d
            * d
            * (
                y[0, 0]
                + y[0, 2]
                + y[2, 0]
                + y[2, 2]
                + 4.0 * (y[0, 1] + y[1, 0] + y[1, 2] + y[2, 1] + 4.0 * y[1, 1])
            )
        )
    if y.shape[0] == 4:  # 1/6, 5/6, 5/6, 1/6
        d = 1.0 / 12.0
        #       return d * ( 5. * ( y[1] + y[2] ) + ( y[0] + y[3] ) )
        w = np.array([1.0, 5.0, 5.0, 1.0])
        ysum = 0.0
        for j in range(0, y.shape[0]):
            for i in range(0, y.shape[1]):
                ysum = ysum + w[i] * w[j] * y[j, i]
        return d * d * ysum
    if y.shape[0] == 5:  # 9/10, 49/90, 64/90, 49/90, 9/90
        d = 1.0 / 180.0
        # return d * ( 64.* y[2] + ( 49. * ( y[1] + y[3] ) )  + 9. * ( y[0] + y[4] ) )
        w = np.array([9.0, 49.0, 64.0, 49.0, 9.0])
        ysum = 0.0
        for j in range(0, y.shape[0]):
            for i in range(0, y.shape[1]):
                ysum = ysum + w[i] * w[j] * y[j, i]
        return d * d * ysum

    raise Exception("Uncoded order")


def lagrange_interp(x, y, q):
    """Lagrange polynomial interpolation. Retruns f(q) which f(x) passes through four data
    points at x[0..3], y[0..3]."""
    # n - numerator, d - denominator
    n0 = (q - x[1]) * (q - x[2]) * (q - x[3])
    d0 = (x[0] - x[1]) * (x[0] - x[2]) * (x[0] - x[3])
    n1 = (q - x[0]) * (q - x[2]) * (q - x[3])
    d1 = (x[1] - x[0]) * (x[1] - x[2]) * (x[1] - x[3])
    n2 = (q - x[0]) * (q - x[1]) * (q - x[3])
    d2 = (x[2] - x[0]) * (x[2] - x[1]) * (x[2] - x[3])
    n3 = (q - x[0]) * (q - x[1]) * (q - x[2])
    d3 = (x[3] - x[0]) * (x[3] - x[1]) * (x[3] - x[2])
    return ((n0 / d0) * y[0] + (n3 / d3) * y[3]) + ((n1 / d1) * y[1] + (n2 / d2) * y[2])


def y_mercator(Ni, phi):
    """Equation (1)"""
    R = Ni / (2 * np.pi)
    return R * (np.log((1.0 + np.sin(phi)) / np.cos(phi)))


def phi_mercator(Ni, y):
    """Equation (2)"""
    R = Ni / (2 * np.pi)
    return np.arctan(np.sinh(y / R)) * (180 / np.pi)  # Converted to degrees


def y_mercator_rounded(Ni, phi):
    y_float = y_mercator(Ni, phi)
    return (np.sign(y_float) * np.round_(np.abs(y_float))).astype(int)


def generate_mercator_grid(
    Ni,
    phi_s,
    phi_n,
    lon0_M,
    lenlon_M,
    refineR,
    shift_equator_to_u_point=True,
    ensure_nj_even=True,
    enhanced_equatorial=False,
):
    print("Requesting Mercator grid with phi range: phi_s,phi_n=", phi_s, phi_n)
    # Diagnose nearest integer y(phi range)
    y_star = y_mercator_rounded(Ni, np.array([phi_s * PI_180, phi_n * PI_180]))
    print("   y*=", y_star, "nj=", y_star[1] - y_star[0] + 1)
    # Ensure that the equator (y=0) is a u-point
    if y_star[0] % 2 == 0:
        print("  *Equator may not be a u-point!")
        # There is another check for this for the whole grid.
        if shift_equator_to_u_point:
            print("  *Fixing this by shifting the bounds!")
            y_star[0] = y_star[0] - 1
            y_star[1] = y_star[1] - 1
            print("   y*=", y_star, "nj=", y_star[1] - y_star[0] + 1)
    if (y_star[1] - y_star[0] + 1) % 2 == 0:
        print("  *Supergrid has an odd number of area cells!")
        if ensure_nj_even:
            print("  *Fixing this by shifting the y_star[1] ")
            y_star[1] = y_star[1] - 1
    Nj = y_star[1] - y_star[0]
    print(
        "   Generating Mercator grid with phi range: phi_s,phi_n=",
        phi_mercator(Ni, y_star),
    )
    phi_M = phi_mercator(Ni, np.arange(y_star[0], y_star[1] + 1))

    # Ensure that the equator (y=0) is included and is a u-point
    equator = 0.0
    equator_index = np.searchsorted(phi_M, equator)
    if equator_index == 0:
        raise Exception("   Ooops: Equator is not in the grid")
    else:
        print("   Equator is at j=", equator_index)
    # Ensure that the equator (y=0) is a u-point
    if equator_index % 2 == 0:
        print("  *Equator is not going to be a u-point of this grid patch.")

    if enhanced_equatorial:
        print("   Enhancing the equator region resolution")
        # Enhance the lattitude resolution between 30S and 30N
        # Set a constant high res lattitude grid spanning 10 degrees centered at the Equator.
        # This construction makes the whole Mercator subgrid symmetric around the Equator.
        #
        # Free MIDAS parameters. Where does this come from and how should it change with resolution?
        phi_enh_d = -5.0  # Starting lattitude of enhanced resolution grid
        phi_cub_d = -30  # Starting lattitude of cubic interpolation

        N_cub = (
            132 * refineR / 2
        )  # Number of points in the cubic interpolation for one shoulder
        # MIDAS has 130, but 132 produces a result closer to 1/2 degree MIDAS grid
        dphi_e = 0.13 * 2 / refineR  # Enhanced resolution 10 degrees around the equator
        N_enh = (
            40 * refineR / 2
        )  # Number of points in the enhanced resolution below equator

        j_c0d = np.where(phi_M < phi_enh_d)[0][
            -1
        ]  # The last index with phi_M<phi_enh_d
        j_phi_cub_d = np.where(phi_M < phi_cub_d)[0][
            -1
        ]  # The last index with phi_M<phi_cub_d
        dphi = phi_M[1:] - phi_M[0:-1]

        cubic_lagrange_interp = True
        cubic_scipy = False

        phi1 = phi_M[0:j_phi_cub_d]
        phi_s = phi_M[j_phi_cub_d - 1]
        dphi_s = phi_M[j_phi_cub_d] - phi_M[j_phi_cub_d - 1]
        phi_e = phi_enh_d

        nodes = [0, 1, N_cub - 2, N_cub - 1]
        phi_nodes = [phi_s, phi_s + dphi_s, phi_e - dphi_e, phi_e]
        q = np.arange(N_cub)

        if cubic_lagrange_interp:
            phi2 = lagrange_interp(nodes, phi_nodes, q)
        elif cubic_scipy:  # MIDAS
            import scipy.interpolate

            f2 = scipy.interpolate.interp1d(nodes, phi_nodes, kind="cubic")
            jInd2 = np.arange(N_cub, dtype=float)
            phi2 = f2(jInd2)

        print(
            "   Meridional range of pure Mercator=(",
            phi1[0],
            ",",
            phi1[-2],
            ") U (",
            -phi1[-2],
            ",",
            -phi1[0],
            ").",
        )
        print(
            "   Meridional range of cubic interpolation=(",
            phi2[0],
            ",",
            phi2[-2],
            ") U (",
            -phi2[-2],
            ",",
            -phi2[0],
            ").",
        )
        phi3 = np.concatenate((phi1[0:-1], phi2))

        phi_s = phi3[-1]
        phi4 = np.linspace(phi_s, 0, int(N_enh))
        print(
            "   Meridional range of enhanced resolution=(", phi4[0], ",", -phi4[0], ")."
        )
        print("   Meridional value of enhanced resolution=", phi4[1] - phi4[0])
        phi5 = np.concatenate((phi3[0:-1], phi4))
        # Make the grid symmetric around the equator!!!!
        phi_M = np.concatenate((phi5[0:-1], -phi5[::-1]))

        # limit the upper lattitude by the requested phi_n
        j_phi_n = np.where(phi_M < phi_n)[0][-1]  # The last index with phi_M<phi_n
        phi_M = phi_M[0:j_phi_n]
        Nj = phi_M.shape[0] - 1

    y_grid_M = np.tile(phi_M.reshape(Nj + 1, 1), (1, Ni + 1))
    lam_M = lon0_M + np.arange(Ni + 1) * lenlon_M / float(Ni)
    x_grid_M = np.tile(lam_M, (Nj + 1, 1))
    # Double check is necessary for enhanced_equatorial
    if y_grid_M.shape[0] % 2 == 0 and ensure_nj_even:
        print(
            "   The number of j's is not even. Fixing this by cutting one row at south."
        )
        y_grid_M = np.delete(y_grid_M, 0, 0)
        x_grid_M = np.delete(x_grid_M, 0, 0)
    print("   Final Mercator grid range=", y_grid_M[0, 0], y_grid_M[-1, 0])
    print("   number of js=", y_grid_M.shape[0])

    return x_grid_M, y_grid_M


###
# Displaced pole cap functions
###
def displacedPoleCap_projection(lon_grid, lat_grid, z_0, r_joint):
    r = np.tan((90 + lat_grid) * PI_180) / r_joint
    # Find the theta that has matching resolution at the unit circle with longitude at the joint
    # This is a conformal transformation of the unit circle (inverse to the one below)
    e2itheta = np.cos(lon_grid * PI_180) + 1j * np.sin(lon_grid * PI_180)
    e2ithetaprime = (e2itheta - z_0) / (1.0 - np.conj(z_0) * e2itheta)
    # Conformal map to displace pole from r=0 to r=r_dispole
    z = r * e2ithetaprime
    w = (z + z_0) / (1 + np.conj(z_0) * z)
    # Inverse projection from tangent plane back to sphere
    lamcDP = np.angle(w, deg=True)
    # lamcDP = np.arctan2(np.imag(w), np.real(w))/PI_180
    # np.angle returns a value in the interval (-180,180)
    # However the input grid longitude is in (-lon0,-lon0+360), e.g., (-300,60)
    # We should shift the angle to be in that interval
    ##But we should also be careful to produce a monotonically increasing longitude, starting from lon0.
    lamcDP = monotonic_bounding(lamcDP, lon_grid[0, 0])
    #
    rw = np.absolute(w)
    phicDP = -90 + np.arctan(rw * r_joint) / PI_180
    return lamcDP, phicDP


def monotonic_bounding(x, x_0):
    x_im1 = x[:, 0] * 0 + x_0  # Initial value
    for i in range(0, x.shape[1]):
        x[:, i] = np.where(x[:, i] - x_im1[:] > 100, x[:, i] - 360, x[:, i])
        x_im1[:] = x[:, i]
    return x


def displacedPoleCap_baseGrid(i, j, ni, nj, lon0, lat0):
    u = lon0 + i * 360.0 / float(ni)
    a = -90.0
    b = lat0
    v = a + j * (b - a) / float(nj)
    du = np.roll(u, shift=-1, axis=0) - u
    dv = np.roll(v, shift=-1, axis=0) - v
    return u, v, du, dv


def displacedPoleCap_mesh(
    i, j, ni, nj, lon0, lat0, lam_pole, r_pole, excluded_fraction=None
):

    long, latg, du, dv = displacedPoleCap_baseGrid(i, j, ni, nj, lon0, lat0)
    lamg = np.tile(long, (latg.shape[0], 1))
    phig = np.tile(latg.reshape((latg.shape[0], 1)), (1, long.shape[0]))
    # Projection from center of globe to plane tangent at south pole
    r_joint = np.tan((90 + lat0) * PI_180)
    z_0 = r_pole * (np.cos(lam_pole * PI_180) + 1j * np.sin(lam_pole * PI_180))
    lams, phis = displacedPoleCap_projection(lamg, phig, z_0, r_joint)

    if excluded_fraction is not None:
        ny, nx = lamg.shape
        jmin = np.ceil(excluded_fraction * ny)
        jmin = jmin + np.mod(jmin, 2)
        jmint = int(jmin)
        return lams[jmint:, :], phis[jmint:, :]
    else:
        return lams, phis


def generate_displaced_pole_grid(Ni, Nj_scap, lon0, lat0, lon_dp, r_dp):
    print("Generating displaced pole grid bounded at latitude ", lat0)
    print("   rdp=", r_dp, " londp=", lon_dp)
    i_s = np.arange(Ni + 1)
    j_s = np.arange(Nj_scap + 1)
    x, y = displacedPoleCap_mesh(i_s, j_s, Ni, Nj_scap, lon0, lat0, lon_dp, r_dp)
    return x, y


# numerical approximation of metrics coefficients h_i and h_j
def great_arc_distance(j0, i0, j1, i1, nx, ny, lon0, lat0, lon_dp, r_dp):
    """Returns great arc distance between nodes (j0,i0) and (j1,i1)"""
    # https://en.wikipedia.org/wiki/Great-circle_distance
    lam0, phi0 = displacedPoleCap_mesh(i0, j0, nx, ny, lon0, lat0, lon_dp, r_dp)
    lam1, phi1 = displacedPoleCap_mesh(i1, j1, nx, ny, lon0, lat0, lon_dp, r_dp)
    lam0, phi0 = lam0 * PI_180, phi0 * PI_180
    lam1, phi1 = lam1 * PI_180, phi1 * PI_180
    dphi, dlam = phi1 - phi0, lam1 - lam0
    # Haversine formula
    d = np.sin(0.5 * dphi) ** 2 + np.sin(0.5 * dlam) ** 2 * np.cos(phi0) * np.cos(phi1)
    return 2.0 * np.arcsin(np.sqrt(d))


def numerical_hi(j, i, nx, ny, lon0, lat0, lon_dp, r_dp, eps, order=6):
    """Returns a numerical approximation to h_lambda"""
    reps = 1.0 / eps
    ds2 = great_arc_distance(j, i + eps, j, i - eps, nx, ny, lon0, lat0, lon_dp, r_dp)
    if order == 2:
        return 0.5 * ds2 * reps
    ds4 = great_arc_distance(
        j, i + 2.0 * eps, j, i - 2.0 * eps, nx, ny, lon0, lat0, lon_dp, r_dp
    )
    if order == 4:
        return (8.0 * ds2 - ds4) * (1.0 / 12.0) * reps
    ds6 = great_arc_distance(
        j, i + 3.0 * eps, j, i - 3.0 * eps, nx, ny, lon0, lat0, lon_dp, r_dp
    )
    if order == 6:
        return (45.0 * ds2 - 9.0 * ds4 + ds6) * (1.0 / 60.0) * reps
    raise Exception("order not coded")


def numerical_hj(j, i, nx, ny, lon0, lat0, lon_dp, r_dp, eps, order=6):
    """Returns a numerical approximation to h_phi"""
    reps = 1.0 / eps
    ds2 = great_arc_distance(j + eps, i, j - eps, i, nx, ny, lon0, lat0, lon_dp, r_dp)
    if order == 2:
        return 0.5 * ds2 * reps
    ds4 = great_arc_distance(
        j + 2.0 * eps, i, j - 2.0 * eps, i, nx, ny, lon0, lat0, lon_dp, r_dp
    )
    if order == 4:
        return (8.0 * ds2 - ds4) * (1.0 / 12.0) * reps
    ds6 = great_arc_distance(
        j + 3.0 * eps, i, j - 3.0 * eps, i, nx, ny, lon0, lat0, lon_dp, r_dp
    )
    if order == 6:
        return (45.0 * ds2 - 9.0 * ds4 + ds6) * (1.0 / 60.0) * reps
    raise Exception("order not coded")


def displacedPoleCap_metrics_quad(
    order, nx, ny, lon0, lat0, lon_dp, r_dp, Re=_default_Re
):
    print("   Calculating displaced pole cap metrics via quadrature ...")
    a, b = quad_positions(order)
    # Note that we need to include the index of the last point of the grid to do the quadrature correctly.
    daq = np.zeros([ny + 1, nx + 1])
    dxq = np.zeros([ny + 1, nx + 1])
    dyq = np.zeros([ny + 1, nx + 1])

    j1d = np.empty([0])
    for j in range(0, ny + 1):
        j_s = b * j + a * (j + 1)
        j1d = np.append(j1d, j_s)

    i1d = np.empty([0])
    for i in range(0, nx + 1):
        i_s = b * i + a * (i + 1)
        i1d = np.append(i1d, i_s)
    # numerical approximation to h_i_in and h_j_inv at quadrature points
    dx = numerical_hi(j1d, i1d, nx, ny, lon0, lat0, lon_dp, r_dp, eps=1e-3, order=order)
    dy = numerical_hj(j1d, i1d, nx, ny, lon0, lat0, lon_dp, r_dp, eps=1e-3, order=order)
    # reshape to send for quad averaging
    dx_r = dx.reshape(ny + 1, order, nx + 1, order)
    dy_r = dy.reshape(ny + 1, order, nx + 1, order)
    # area element
    dxdy_r = dx_r * dy_r

    for j in range(0, ny + 1):
        for i in range(0, nx + 1):
            daq[j, i] = quad_average_2d(dxdy_r[j, :, i, :])
            dxq[j, i] = quad_average(dx_r[j, 0, i, :])
            dyq[j, i] = quad_average(dy_r[j, :, i, 0])

    daq = daq[:-1, :-1] * Re * Re
    dxq = dxq[:, :-1] * Re
    dyq = dyq[:-1, :] * Re

    return dxq, dyq, daq


def cut_below(lam, phi, lowerlat):
    nj, ni = lam.shape
    for j in range(0, nj):
        if phi[j, 0] > lowerlat:
            break
    jmin = j
    #    print("jmin",jmin)
    return lam[jmin:, :], phi[jmin:, :]


def cut_above(lam, phi, upperlat):
    nj, ni = lam.shape
    for j in range(0, nj):
        if phi[j, 0] > upperlat:
            break
    jmax = j
    #    print("jmax",jmax)
    return lam[0:jmax, :], phi[0:jmax, :]


# utility function to plot grids
def plot_mesh_in_latlon(
    lam, phi, stride=1, phi_color="k", lam_color="r", newfig=True, title=None
):
    import matplotlib.pyplot as plt

    #    import seaborn as sns; sns.set()
    if phi.shape != lam.shape:
        raise Exception("Ooops: lam and phi should have same shape")
    nj, ni = lam.shape
    if newfig:
        plt.figure(figsize=(10, 10))
    for i in range(0, ni, stride):
        plt.plot(lam[:, i], phi[:, i], lam_color)
    for j in range(0, nj, stride):
        plt.plot(lam[j, :], phi[j, :], phi_color)
    if title is not None:
        plt.title(title)


#    plt.show()


def plot_mesh_in_xyz(
    lam,
    phi,
    stride=1,
    phi_color="k",
    lam_color="r",
    lowerlat=None,
    upperlat=None,
    newfig=True,
    title=None,
):
    if lowerlat is not None:
        lam, phi = cut_below(lam, phi, lowerlat=lowerlat)
    if upperlat is not None:
        lam, phi = cut_above(lam, phi, upperlat=upperlat)
    x = np.cos(phi * PI_180) * np.cos(lam * PI_180)
    y = np.cos(phi * PI_180) * np.sin(lam * PI_180)
    z = np.sin(phi * PI_180)
    plot_mesh_in_latlon(
        x,
        y,
        stride=stride,
        phi_color=phi_color,
        lam_color=lam_color,
        newfig=newfig,
        title=title,
    )


def displacedPoleCap_plot(x_s, y_s, lon0, lon_dp, lat0_SO, stride=40, block=False):
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 10))
    x_s = x_s + lon0
    plot_mesh_in_xyz(x_s, y_s, stride=stride, newfig=False)
    lam__ = x_s[-1, 0]
    phi__ = y_s[-1, 0]
    x__ = np.cos(phi__ * PI_180) * np.cos(lam__ * PI_180)
    y__ = np.cos(phi__ * PI_180) * np.sin(lam__ * PI_180)
    plt.annotate(
        "i=0",
        xy=(x__, y__),
        xytext=(x__ - 0.01, y__ + 0.01),
        arrowprops=dict(arrowstyle="->"),
    )

    # plot the location of the displaced pole
    lam__ = np.array([lon_dp, lon_dp]) + lon0
    phi__ = np.array([lat0_SO + 2, lat0_SO - 2])
    x__ = np.cos(phi__ * PI_180) * np.cos(lam__ * PI_180)
    y__ = np.cos(phi__ * PI_180) * np.sin(lam__ * PI_180)
    plt.plot(x__, y__, linestyle="dashed")
    plt.annotate(
        "pole's lon",
        xy=(x__[0], y__[0]),
        xytext=(x__[0], y__[0] + 0.01),
        arrowprops=dict(arrowstyle="->"),
    )

    n_i = x_s.shape[1] - 1
    poles_i = int(n_i * np.mod(lon_dp - lon0, 360) / 360.0)
    lam__ = x_s[:, poles_i]
    phi__ = y_s[:, poles_i]
    x__ = np.cos(phi__ * PI_180) * np.cos(lam__ * PI_180)
    y__ = np.cos(phi__ * PI_180) * np.sin(lam__ * PI_180)
    plt.plot(x__, y__, linewidth=4, color="green")

    antipoles_i = int(n_i * np.mod(180 + lon_dp - lon0, 360) / 360.0)
    lam__ = x_s[:, antipoles_i]
    phi__ = y_s[:, antipoles_i]
    x__ = np.cos(phi__ * PI_180) * np.cos(lam__ * PI_180)
    y__ = np.cos(phi__ * PI_180) * np.sin(lam__ * PI_180)
    plt.plot(x__, y__, linewidth=4, color="blue")
    plt.show(block)


def displacedPoleCap_plot2(x_s, y_s, lon0, lon_dp, lat0, stride=40, block=False):
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 10))
    m = Basemap(projection="spstere", boundinglat=-60, lon_0=300, resolution="l")
    m.drawcoastlines()
    m.fillcontinents(color="coral", lake_color="aqua")
    m.drawparallels(np.arange(lat0, 0.0, 10.0), labels=[True, True])
    m.drawmeridians(np.arange(0.0, 360.0, 20.0), labels=[True, True, True, True])
    m.drawmeridians(
        np.arange(lon_dp, lon_dp + 1, 1.0),
        labels=[True, True, True, True],
        color="r",
        linewidth=4,
    )
    m.drawmapboundary(fill_color="aqua")

    x, y = m(x_s[::2, ::80], y_s[::2, ::80])
    plt.plot(x, y, "ok", markersize=2)
    plt.show(block)


def mdist(x1, x2):
    """Returns positive distance modulo 360."""
    return np.minimum(np.mod(x1 - x2, 360.0), np.mod(x2 - x1, 360.0))


def generate_grid_metrics_MIDAS(
    x, y, axis_units="degrees", Re=_default_Re, latlon_areafix=True
):
    nytot, nxtot = x.shape
    if axis_units == "m":
        metric = 1.0
    if axis_units == "km":
        metric = 1.0e3
    if axis_units == "degrees":
        metric = Re * PI_180
    lv = (0.5 * (y[:, 1:] + y[:, :-1])) * PI_180
    dx_i = mdist(x[:, 1:], x[:, :-1]) * PI_180
    dy_i = (y[:, 1:] - y[:, :-1]) * PI_180
    dx = Re * np.sqrt(dy_i ** 2 + (dx_i * np.cos(lv)) ** 2)
    lu = (0.5 * (y[1:, :] + y[:-1, :])) * PI_180
    dx_j = mdist(x[1:, :], x[:-1, :]) * PI_180
    dy_j = (y[1:, :] - y[:-1, :]) * PI_180
    dy = Re * np.sqrt(dy_j ** 2 + (dx_j * np.cos(lu)) ** 2)

    ymid_j = 0.5 * (y + np.roll(y, shift=-1, axis=0))
    ymid_i = 0.5 * (y + np.roll(y, shift=-1, axis=1))
    dy_j = np.roll(y, shift=-1, axis=0) - y
    dy_i = np.roll(y, shift=-1, axis=1) - y
    dx_i = mdist(np.roll(x, shift=-1, axis=1), x)
    dx_j = mdist(np.roll(x, shift=-1, axis=0), x)
    if latlon_areafix:
        sl = np.sin(lv)
        dx_i = mdist(x[:, 1:], x[:, :-1]) * PI_180
        area = (Re ** 2) * (
            (0.5 * (dx_i[1:, :] + dx_i[:-1, :])) * (sl[1:, :] - sl[:-1, :])
        )
    else:
        area = 0.25 * ((dx[1:, :] + dx[:-1, :]) * (dy[:, 1:] + dy[:, :-1]))
    return dx, dy, area


def angle_x(x, y):
    """Returns the orientation angle of the grid box"""
    if x.shape != y.shape:
        raise Exception("Input arrays do not have the same shape!")
    angle_dx = np.zeros(x.shape)
    # The corrected version of angle_dx, in addition to including spherical metrics, is centered in the interior and one-sided at the grid edges
    angle_dx[:, 1:-1] = np.arctan2(
        y[:, 2:] - y[:, :-2], (x[:, 2:] - x[:, :-2]) * np.cos(y[:, 1:-1] * PI_180)
    )
    angle_dx[:, 0] = np.arctan2(
        y[:, 1] - y[:, 0], (x[:, 1] - x[:, 0]) * np.cos(y[:, 0] * PI_180)
    )
    angle_dx[:, -1] = np.arctan2(
        y[:, -1] - y[:, -2], (x[:, -1] - x[:, -2]) * np.cos(y[:, -1] * PI_180)
    )
    angle_dx = angle_dx / PI_180
    return angle_dx


def metrics_error(
    dx_,
    dy_,
    area_,
    Ni,
    lat1,
    lat2=90,
    Re=_default_Re,
    bipolar=False,
    displaced_pole=-999,
    excluded_fraction=None,
):
    exact_area = (
        2 * np.pi * (Re ** 2) * np.abs(np.sin(lat2 * PI_180) - np.sin(lat1 * PI_180))
    )
    exact_lat_arc_length = np.abs(lat2 - lat1) * PI_180 * Re
    exact_lon_arc_length = np.cos(lat1 * PI_180) * 2 * np.pi * Re
    grid_lat_arc_length = np.sum(dy_[:, Ni // 4])
    grid_lon_arc_length = np.sum(dx_[0, :])
    if lat1 > lat2:
        grid_lon_arc_length = np.sum(dx_[-1, :])
    if bipolar:
        # length of the fold
        grid_lon_arc_length2 = np.sum(dx_[-1, :])
        # This must be 4*grid_lat_arc_length
        lon_arc2_error = (
            100
            * (grid_lon_arc_length2 / 4 - exact_lat_arc_length)
            / exact_lat_arc_length
        )
    area_error = 100 * (np.sum(area_) - exact_area) / exact_area
    lat_arc_error = (
        100 * (grid_lat_arc_length - exact_lat_arc_length) / exact_lat_arc_length
    )
    lon_arc_error = (
        100 * (grid_lon_arc_length - exact_lon_arc_length) / exact_lon_arc_length
    )
    if displaced_pole != -999:
        antipole = displaced_pole + Ni // 2
        if displaced_pole > Ni // 2:
            antipole = displaced_pole - Ni // 2
        grid_lat_arc_length = np.sum(dy_[:, displaced_pole]) + np.sum(dy_[:, antipole])
        lat_arc_error = (
            100
            * (grid_lat_arc_length - 2.0 * exact_lat_arc_length)
            / exact_lat_arc_length
        )
    if excluded_fraction:
        print(
            "   Cannot estimate area and dy accuracies with excluded_fraction (doughnut)! "
        )
    if bipolar:
        return area_error, lat_arc_error, lon_arc_error, lon_arc2_error
    else:
        return area_error, lat_arc_error, lon_arc_error


def write_nc(
    x,
    y,
    dx,
    dy,
    area,
    angle_dx,
    axis_units="degrees",
    fnam=None,
    format="NETCDF3_64BIT",
    description=None,
    history=None,
    source=None,
    no_changing_meta=None,
    debug=False,
):
    import netCDF4 as nc

    if fnam is None:
        fnam = "supergrid.nc"
    fout = nc.Dataset(fnam, "w", clobber=True, format=format)

    if debug:
        chksum(x, "x")
        chksum(y, "y")
        chksum(dx, "dx")
        chksum(dy, "dy")
        chksum(area, "area")
        chksum(angle_dx, "angle_dx")

    ny = area.shape[0]
    nx = area.shape[1]
    nyp = ny + 1
    nxp = nx + 1
    print("   Writing netcdf file with ny,nx= ", ny, nx)

    nyp = fout.createDimension("nyp", nyp)
    nxp = fout.createDimension("nxp", nxp)
    ny = fout.createDimension("ny", ny)
    nx = fout.createDimension("nx", nx)
    string = fout.createDimension("string", 255)
    tile = fout.createVariable("tile", "S1", ("string"))
    yv = fout.createVariable("y", "f8", ("nyp", "nxp"))
    xv = fout.createVariable("x", "f8", ("nyp", "nxp"))
    yv.units = "degrees"
    xv.units = "degrees"
    yv[:] = y
    xv[:] = x
    stringvals = np.empty(1, "S" + repr(len(tile)))
    stringvals[0] = "tile1"
    tile[:] = nc.stringtochar(stringvals)
    dyv = fout.createVariable("dy", "f8", ("ny", "nxp"))
    dyv.units = "meters"
    dyv[:] = dy
    dxv = fout.createVariable("dx", "f8", ("nyp", "nx"))
    dxv.units = "meters"
    dxv[:] = dx
    areav = fout.createVariable("area", "f8", ("ny", "nx"))
    areav.units = "m2"
    areav[:] = area
    anglev = fout.createVariable("angle_dx", "f8", ("nyp", "nxp"))
    anglev.units = "degrees"
    anglev[:] = angle_dx
    # global attributes
    if not no_changing_meta:
        fout.history = history
        fout.description = description
        fout.source = source

    fout.sync()
    fout.close()


def generate_latlon_grid(
    lni, lnj, llon0, llen_lon, llat0, llen_lat, ensure_nj_even=True
):
    print("Generating regular lat-lon grid between latitudes ", llat0, llat0 + llen_lat)
    llonSP = llon0 + np.arange(lni + 1) * llen_lon / float(lni)
    llatSP = llat0 + np.arange(lnj + 1) * llen_lat / float(lnj)
    if llatSP.shape[0] % 2 == 0 and ensure_nj_even:
        print(
            "   The number of j's is not even. Fixing this by cutting one row at south."
        )
        llatSP = np.delete(llatSP, 0, 0)

    llamSP = np.tile(llonSP, (llatSP.shape[0], 1))
    lphiSP = np.tile(llatSP.reshape((llatSP.shape[0], 1)), (1, llonSP.shape[0]))

    print(
        "   generated regular lat-lon grid between latitudes ",
        lphiSP[0, 0],
        lphiSP[-1, 0],
    )
    print("   number of js=", lphiSP.shape[0])

    #    h_i_inv=llen_lon*PI_180*np.cos(lphiSP*PI_180)/lni
    #    h_j_inv=llen_lat*PI_180*np.ones(lphiSP.shape)/lnj
    #    delsin_j = np.roll(np.sin(lphiSP*PI_180),shift=-1,axis=0) - np.sin(lphiSP*PI_180)
    #    dx_h=h_i_inv[:,:-1]*_default_Re
    #    dy_h=h_j_inv[:-1,:]*_default_Re
    #    area=delsin_j[:-1,:-1]*_default_Re*_default_Re*llen_lon*PI_180/lni

    return llamSP, lphiSP


def usage():
    print(
        "ocean_grid_generator.py -f <output_grid_filename> -r <inverse_degrees_resolution> [--rdp=<displacement_factor/0.2> --south_cutoff_ang=<degrees_south_to_start> --south_cutoff_row=<rows_south_to_cut> --reproduce_MIDAS_grids --match_dy --even_j --plot --write_subgrid_files --enhanced_equatorial --no-metrics --gridlist=sc]"
    )


def main(argv):

    degree_resolution_inverse = (
        4  # (2 for half) or (4 for quarter) or (8 for 1/8) degree grid
    )
    south_cap = True
    gridfilename = "tripolar_res" + str(degree_resolution_inverse) + ".nc"
    r_dp = 0.0  # r value   of the displaced pole
    lon_dp = 80.0  # longitude of the displaced pole
    south_cutoff_row = 0
    south_cutoff_ang = -90.0
    reproduce_MIDAS_grids = False
    match_dy = False
    write_subgrid_files = False
    plotem = False
    no_changing_meta = False
    enhanced_equatorial = False
    debug = False
    grids = "bipolar,mercator,so,sc,all"
    calculate_metrics = True
    # Ensure the number of j partitions are even for the sub-grids
    ensure_nj_even = False
    shift_equator_to_u_point = True

    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "hdf:r:",
            [
                "gridfilename=",
                "inverse_resolution=",
                "south_cutoff_ang=",
                "south_cutoff_row=",
                "rdp=",
                "londp=",
                "reproduce_MIDAS_grids",
                "match_dy",
                "even_j",
                "plot",
                "write_subgrid_files",
                "no_changing_meta",
                "enhanced_equatorial",
                "no-metrics",
                "gridlist=",
            ],
        )
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == "-h":
            usage()
            sys.exit()
        elif opt in ("-f", "--gridfilename"):
            gridfilename = arg
        elif opt in ("-r", "--inverse_resolution"):
            degree_resolution_inverse = float(arg)
        elif opt in ("-d"):
            debug = True
        elif opt in ("--south_cutoff_ang"):
            south_cutoff_ang = float(arg)
        elif opt in ("--south_cutoff_row"):
            south_cutoff_row = int(arg)
        elif opt in ("--rdp"):
            r_dp = float(arg)
        elif opt in ("--londp"):
            lon_dp = float(arg)
        elif opt in ("--reproduce_MIDAS_grids"):
            reproduce_MIDAS_grids = True
        elif opt in ("--match_dy"):
            match_dy = True
        elif opt in ("--even_j"):
            ensure_nj_even = True
        elif opt in ("--plot"):
            plotem = True
        elif opt in ("--no-metrics"):
            calculate_metrics = False
        elif opt in ("--write_subgrid_files"):
            write_subgrid_files = True
        elif opt in ("--no_changing_meta"):
            no_changing_meta = True
        elif opt in ("--enhanced_equatorial"):
            enhanced_equatorial = True
        elif opt in ("--gridlist"):
            grids = arg
        else:
            assert False, "unhandled option"

    # Information to write in file as metadata
    if not no_changing_meta:
        import socket

        host = str(socket.gethostname())
        scriptpath = sys.argv[0]
        scriptbasename = (
            subprocess.check_output("basename " + scriptpath, shell=True)
            .decode("ascii")
            .rstrip("\n")
        )
        scriptdirname = (
            subprocess.check_output("dirname " + scriptpath, shell=True)
            .decode("ascii")
            .rstrip("\n")
        )
        scriptgithash = (
            subprocess.check_output(
                "cd " + scriptdirname + ";git rev-parse HEAD; exit 0",
                stderr=subprocess.STDOUT,
                shell=True,
            )
            .decode("ascii")
            .rstrip("\n")
        )
        scriptgitMod = (
            subprocess.check_output(
                "cd "
                + scriptdirname
                + ";git status --porcelain "
                + scriptbasename
                + " | awk '{print $1}' ; exit 0",
                stderr=subprocess.STDOUT,
                shell=True,
            )
            .decode("ascii")
            .rstrip("\n")
        )
        if "M" in str(scriptgitMod):
            scriptgitMod = " , But was localy Modified!"

    hist = "This grid file was generated via command " + " ".join(sys.argv)
    if not no_changing_meta:
        hist = hist + " on " + str(datetime.date.today()) + " on platform " + host

    desc = (
        "This is an orthogonal coordinate grid for the Earth with a nominal resoution of "
        + str(1 / degree_resolution_inverse)
        + " degrees along the equator. "
    )

    source = ""
    if not no_changing_meta:
        source = source + scriptpath + " had git hash " + scriptgithash + scriptgitMod
        source = (
            source
            + ". To obtain the grid generating code do: git clone  https://github.com/nikizadehgfdl/grid_generation.git ; cd grid_generation;  git checkout "
            + scriptgithash
        )

    # Specify the default grid properties
    refineS = 2  # factor 2 is for supergrid
    refineR = degree_resolution_inverse
    lenlon = 360  # global longitude range
    lon0 = -300.0  # Starting longitude of the map
    Ni = int(refineR * refineS * lenlon)
    ###
    ###Mercator grid
    ###
    # MIDAS has nominal starting latitude for Mercator grid = -65 for 1/4 degree, -70 for 1/2 degree
    # MIDAS has nominal latitude range of Mercator grid     = 125 for 1/4 degree, 135 for 1/2 degree
    # Instead we use:
    phi_s_Merc, phi_n_Merc = -66.85954725, 64.05895973
    if refineR == 2:
        # phi_s_Merc, phi_n_Merc = -68.05725376601046, 65.0 #These give a 1/2 degree enhanced equatorial close to MIDAS result
        # shift_equator_to_u_point= False
        phi_s_Merc, phi_n_Merc = -68.0, 65.0
    ###
    # Southern Ocean grid
    ###
    lat0_SO = -78.0  # Starting lat of Southern Ocean grid
    lenlat_SO = phi_s_Merc - lat0_SO
    deltaPhiSO = 1.0 / refineR / refineS
    # To get the same number of points as existing 1/2 and 1/4 degree grids that were generated with MIDAS
    Nj_SO = int(refineR * 55)
    if refineR == 2:
        Nj_SO = 54 * refineS + 1
    ###
    # Bipolar cap
    ###
    lon_bp = lon0  # longitude of the bipole(s)
    # To get the same number of points as existing 1/2 and 1/4 degree grids that were generated with MIDAS
    Nj_ncap = int(60 * refineR * refineS)
    if refineR == 2:
        Nj_ncap = 119 * refineS
    ###
    # South cap
    ###
    # To get the same number of points as existing 1/4 degree grids that were generated with MIDAS
    Nj_scap = int(refineR * 40)
    # Refine the grid by a factor and then exclude the inner circle corresponding to that factor
    # These factors are heuristic and we adjust them to get a grid with the same number of points
    # as the existing 1/4 degree grid of OM4p25
    Nj_scap = Nj_scap * 7 // 4
    doughnut = 0.28 * 7 / 4
    lat0_SC = lat0_SO

    if "mercator" in grids or "all" in grids:
        if not reproduce_MIDAS_grids:
            lamMerc, phiMerc = generate_mercator_grid(
                Ni,
                phi_s_Merc,
                phi_n_Merc,
                lon0,
                lenlon,
                refineR,
                shift_equator_to_u_point=shift_equator_to_u_point,
                ensure_nj_even=ensure_nj_even,
                enhanced_equatorial=enhanced_equatorial,
            )
            angleMerc = angle_x(lamMerc, phiMerc)
            dxMerc = -np.ones([lamMerc.shape[0], lamMerc.shape[1] - 1])
            dyMerc = -np.ones([lamMerc.shape[0] - 1, lamMerc.shape[1]])
            areaMerc = -np.ones([lamMerc.shape[0] - 1, lamMerc.shape[1] - 1])
            if calculate_metrics:
                # For spherical grids we can safely use the MIDAS algorithm for calculating the metrics
                dxMerc, dyMerc, areaMerc = generate_grid_metrics_MIDAS(lamMerc, phiMerc)
                print(
                    "   CHECK_metrics_MIDAS: % errors in (area, lat arc, lon arc)",
                    metrics_error(
                        dxMerc, dyMerc, areaMerc, Ni, phiMerc[0, 0], phiMerc[-1, 0]
                    ),
                )
        else:  # use pymidas package
            from pymidas.rectgrid_gen import supergrid

            if refineR == 2:
                Nj_Merc = 364 * refineS
                phi_s_Merc = (
                    -70.0
                )  # This is a nominal starting latitude for Mercator grid
                lenlat_Merc = 135.0  # nominal latitude range of Mercator grid
            elif refineR == 4:
                Nj_Merc = 700 * refineS
                phi_s_Merc = (
                    -65.0
                )  # This is a nominal starting latitude for Mercator grid
                lenlat_Merc = 125.0  # nominal latitude range of Mercator grid
            else:
                raise Exception("Unknown resolution for --reproduce_MIDAS_grids")
            #### Begin Mercator Grid
            print("constructing a mercator supergrid with (ni,nj) = ", Ni, Nj_Merc)
            print("nominal starting lat and starting longitude =", phi_s_Merc, lon0)
            print("and nominal width in latitude = ", lenlat_Merc)
            mercator = supergrid(
                Ni,
                Nj_Merc,
                "mercator",
                "degrees",
                phi_s_Merc,
                lenlat_Merc,
                lon0,
                360.0,
                cyclic_x=True,
            )

            # Add equatorial enhancement for 1/2 degree MIDAS grid
            if refineR == 2:
                print("Enhancing the equator resolution")
                import scipy.interpolate

                phi = np.ascontiguousarray(mercator.y[:, 0])
                dphi = phi[1:] - phi[0:-1]
                phi = mercator.y[:, 0]
                dphi = phi[1:] - phi[0:-1]
                jind = np.where(phi > -30.0)[0][0]
                jind = jind + np.mod(jind, 2)
                phi = 1.0 * phi[0:jind]
                dphi = dphi[0:jind]
                N = 130
                phi_s = phi[-1]
                dphi_s = dphi[-1]
                phi_e = -5.0
                dphi_e = 0.13
                nodes = [0, 1, N - 2, N - 1]
                phi_nodes = [phi_s, phi_s + dphi_s, phi_e - dphi_e, phi_e]
                f2 = scipy.interpolate.interp1d(nodes, phi_nodes, kind="cubic")
                jInd2 = np.arange(N, dtype=float)
                phi2 = f2(jInd2)
                print(
                    "Meridional range of pure Mercator=(",
                    phi[0],
                    ",",
                    phi[-2],
                    ") U (",
                    -phi[-2],
                    ",",
                    -phi[0],
                    ").",
                )
                print(
                    "Meridional range of cubic interpolation=(",
                    phi2[0],
                    ",",
                    phi2[-2],
                    ") U (",
                    -phi2[-2],
                    ",",
                    -phi2[0],
                    ").",
                )
                phi = np.concatenate((phi[0:-1], phi2))
                N = 40
                phi_s = phi[-1]
                phi2 = np.linspace(phi_s, 0, N)
                print(
                    "Meridional range of enhanced resolution=(",
                    phi2[0],
                    ",",
                    -phi2[0],
                    ").",
                )
                print("Meridional value of enhanced resolution=", phi2[1] - phi2[0])
                PHI = np.concatenate((phi[0:-1], phi2))
                PHI = np.concatenate((PHI[0:-1], -PHI[::-1]))
                LAMBDA = np.linspace(lon0, lon0 + 360.0, Ni + 1)
                jind = np.where(PHI > -78.0)[0][0]
                jind = jind + np.mod(jind, 2)
                jind2 = np.where(PHI > 65.0)[0][0]
                jind2 = jind2 + np.mod(jind2, 2)
                PHI2 = PHI[jind : jind2 - 1]
                x, y = np.meshgrid(LAMBDA, PHI2)
                mercator = supergrid(
                    xdat=x, ydat=y, axis_units="degrees", cyclic_x=True
                )

            print("mercator.y.shape= ", mercator.y.shape)
            print("mercator max/min latitude=", mercator.y.max(), mercator.y.min())
            print("mercator start/end longitude=", mercator.x[0, 0], mercator.x[0, -1])
            print("mercator start/end latitudes=", mercator.y[0, 0], mercator.y[-1, 0])
            mercator.grid_metrics()
            lamMerc, phiMerc = mercator.x, mercator.y
            dxMerc, dyMerc, areaMerc, angleMerc = (
                mercator.dx,
                mercator.dy,
                mercator.area,
                mercator.angle_dx,
            )
            print(
                "   CHECK_metrics_MIDAS: % errors in (area, lat arc, lon arc)",
                metrics_error(
                    dxMerc, dyMerc, areaMerc, Ni, phiMerc[0, 0], phiMerc[-1, 0]
                ),
            )

        if write_subgrid_files:
            write_nc(
                lamMerc,
                phiMerc,
                dxMerc,
                dyMerc,
                areaMerc,
                angleMerc,
                axis_units="degrees",
                fnam=gridfilename + "Merc.nc",
                description=desc,
                history=hist,
                source=source,
                debug=debug,
            )

        # The phi resolution in the first and last row of Mercator grid along the symmetry meridian
        DeltaPhiMerc_so = phiMerc[1, Ni // 4] - phiMerc[0, Ni // 4]
        DeltaPhiMerc_no = phiMerc[-1, Ni // 4] - phiMerc[-2, Ni // 4]
        # Start lattitude from dy above the last Mercator grid
        lat0_bp = phiMerc[-1, Ni // 4] + DeltaPhiMerc_no

        if match_dy:
            # The bopolar grid should start from the same lattitude that Mercator ends.
            # Then when we combine the two grids we should drop the x,y,dx,angle from one of the two.
            # This way we get a continous dy and area.
            lat0_bp = phiMerc[-1, Ni // 4]
            # Determine the number of bipolar cap grid point in the y direction such that the y resolution
            # along symmetry meridian is a constant and is equal to (continuous with) the last Mercator dy.
            # Note that int(0.5+x) is used to return the nearest integer to a float with deterministic
            # behavior for middle points.
            # Note that int(0.5+x) is equivalent to math.floor(0.5+x)
            Nj_ncap = int(
                0.5 + (90.0 - lat0_bp) / DeltaPhiMerc_no
            )  # Impose boundary condition for smooth dy

            # Make the last SO grid point a (Mercator) step below the first Mercator lattitude.
            # lenlat_SO = phiMerc[0,Ni//4] - DeltaPhiMerc_so - lat0_SO #Start from a lattitude to smooth out dy.
            # Niki: I think this is wrong!
            # The SO grid should end at the same lattitude that Mercator starts.
            # Then when we combine the two grids we should drop the x,y,dx,angle from one of the two.
            # This way we get a continous dy and area.
            lenlat_SO = phiMerc[0, Ni // 4] - lat0_SO
            # Determine the number of grid point in the y direction such that the y resolution is equal to
            # (continuous with) the first Mercator dy.
            Nj_SO = int(
                0.5 + lenlat_SO / DeltaPhiMerc_so
            )  # Make the resolution continious with the Mercator at joint

    ###
    ###Northern bipolar cap
    ###
    if "bipolar" in grids or "all" in grids:
        if not reproduce_MIDAS_grids:
            # Generate the bipolar grid
            lamBP, phiBP, dxBP_h, dyBP_h = generate_bipolar_cap_mesh(
                Ni, Nj_ncap, lat0_bp, lon_bp, ensure_nj_even=ensure_nj_even
            )
            # Metrics via quadratue of h's
            rp = np.tan(0.5 * (90 - lat0_bp) * PI_180)
            dxBP = -np.ones([lamBP.shape[0], lamBP.shape[1] - 1])
            dyBP = -np.ones([lamBP.shape[0] - 1, lamBP.shape[1]])
            areaBP = -np.ones([lamBP.shape[0] - 1, lamBP.shape[1] - 1])
            if calculate_metrics:
                dxBP, dyBP, areaBP = bipolar_cap_metrics_quad_fast(
                    5, phiBP.shape[1] - 1, phiBP.shape[0] - 1, lat0_bp, lon_bp, rp
                )
                print(
                    "   CHECK_metrics_hquad: % errors in (area, lat arc, lon arc1, lon arc2)",
                    metrics_error(dxBP, dyBP, areaBP, Ni, lat0_bp, 90.0, bipolar=True),
                )
            angleBP = angle_x(lamBP, phiBP)
        else:
            lat0_bp = mercator.y.max()
            dlat = 90.0 - lat0_bp

            tripolar_n = supergrid(
                Ni,
                Nj_ncap,
                "spherical",
                "degrees",
                lat0_bp,
                dlat,
                lon0,
                360.0,
                tripolar_n=True,
            )
            tripolar_n.grid_metrics()
            print(
                "generated a tripolar supergrid of size (ny,nx)= ",
                tripolar_n.y.shape[0] - 1,
                tripolar_n.y.shape[1] - 1,
            )
            print("tripolar grid starting longitude = ", tripolar_n.x[0, 0])
            print("tripolar grid starting latitude = ", tripolar_n.y[0, 0])

            lamBP, phiBP = tripolar_n.x, tripolar_n.y
            dxBP, dyBP, areaBP, angleBP = (
                tripolar_n.dx,
                tripolar_n.dy,
                tripolar_n.area,
                tripolar_n.angle_dx,
            )

            print(
                "   CHECK_metrics_MIDAS: % errors in (area, lat arc, lon arc)",
                metrics_error(dxBP, dyBP, areaBP, Ni, lat0_bp, 90.0),
            )

        if write_subgrid_files:
            write_nc(
                lamBP,
                phiBP,
                dxBP,
                dyBP,
                areaBP,
                angleBP,
                axis_units="degrees",
                fnam=gridfilename + "BP.nc",
                description=desc,
                history=hist,
                source=source,
                debug=debug,
            )

    if "so" in grids or "all" in grids:
        ###
        ###Southern Ocean grid
        ###
        if not reproduce_MIDAS_grids:
            lamSO, phiSO = generate_latlon_grid(
                Ni,
                Nj_SO,
                lon0,
                lenlon,
                lat0_SO,
                lenlat_SO,
                ensure_nj_even=ensure_nj_even,
            )
            dxSO = -np.ones([lamSO.shape[0], lamSO.shape[1] - 1])
            dySO = -np.ones([lamSO.shape[0] - 1, lamSO.shape[1]])
            areaSO = -np.ones([lamSO.shape[0] - 1, lamSO.shape[1] - 1])
            if calculate_metrics:
                # For spherical grids we can safely use the MIDAS algorithm for calculating the metrics
                dxSO, dySO, areaSO = generate_grid_metrics_MIDAS(lamSO, phiSO)
            angleSO = angle_x(lamSO, phiSO)
            print(
                "   CHECK_metrics_MIDAS: % errors in (area, lat arc, lon arc)",
                metrics_error(dxSO, dySO, areaSO, Ni, phiSO[0, 0], phiSO[-1, 0]),
            )
        else:
            print("constructing a spherical supergrid with (ny,nx) = ", Ni, Nj_SO)
            print("nominal starting lat and starting longitude =", lat0_SO, lon0)
            print("and nominal width in latitude = ", mercator.y.min() - lat0_SO)
            spherical = supergrid(
                Ni,
                Nj_SO,
                "spherical",
                "degrees",
                lat0_SO,
                mercator.y.min() - lat0_SO,
                lon0,
                360.0,
                cyclic_x=True,
            )
            spherical.grid_metrics()
            print(
                "southern ocean spherical max/min latitude=",
                spherical.y.max(),
                spherical.y.min(),
            )
            print(
                "spherical nj,ni=", spherical.y.shape[0] - 1, spherical.y.shape[1] - 1
            )
            print("spherical starting longitude=", spherical.x[0, 0])
            print("spherical ending longitude=", spherical.x[0, -1])

            lamSO, phiSO = spherical.x, spherical.y
            dxSO, dySO, areaSO, angleSO = (
                spherical.dx,
                spherical.dy,
                spherical.area,
                spherical.angle_dx,
            )
            # Metrics errors via MIDAS
            print(
                "   CHECK_metrics_MIDAS: % errors in (area, lat arc, lon arc)",
                metrics_error(dxSO, dySO, areaSO, Ni, phiSO[0, 0], phiSO[-1, 0]),
            )

        deltaPhiSO = phiSO[1, Ni // 4] - phiSO[0, Ni // 4]
        lat0_SC = phiSO[0, Ni // 4] - deltaPhiSO
        # The above heuristics produce a displaced pole grid with a nominal resolution
        # different from the rest of the grid!
        # To get the nominal resolution right  we must instead make the resolution continuous across the joint
        # fullArc = lat0_SC+90.
        # if(match_dy):
        #    Nj_scap = int(fullArc/deltaPhiSO)

        if write_subgrid_files:
            write_nc(
                lamSO,
                phiSO,
                dxSO,
                dySO,
                areaSO,
                angleSO,
                axis_units="degrees",
                fnam=gridfilename + "SO.nc",
                description=desc,
                history=hist,
                source=source,
                debug=debug,
            )

    if "sc" in grids or "all" in grids:
        ###
        ###Southern cap
        ###
        if r_dp == 0.0:
            fullArc = lat0_SC + 90.0
            Nj_scap = int(fullArc / deltaPhiSO)
            if not reproduce_MIDAS_grids:
                lamSC, phiSC = generate_latlon_grid(
                    Ni,
                    Nj_scap,
                    lon0,
                    lenlon,
                    -90.0,
                    90 + lat0_SO,
                    ensure_nj_even=ensure_nj_even,
                )
                angleSC = angle_x(lamSC, phiSC)
                dxSC = -np.ones([lamSC.shape[0], lamSC.shape[1] - 1])
                dySC = -np.ones([lamSC.shape[0] - 1, lamSC.shape[1]])
                areaSC = -np.ones([lamSC.shape[0] - 1, lamSC.shape[1] - 1])
                if calculate_metrics:
                    # For spherical grids we can safely use the MIDAS algorithm for calculating the metrics
                    dxSC, dySC, areaSC = generate_grid_metrics_MIDAS(lamSC, phiSC)
                    print(
                        "   CHECK_metrics_MIDAS: % errors in (area, lat arc, lon arc)",
                        metrics_error(
                            dxSC, dySC, areaSC, Ni, phiSC[-1, 0], phiSC[0, 0]
                        ),
                    )
            else:
                spherical_cap = supergrid(
                    Ni,
                    Nj_scap,
                    "spherical",
                    "degrees",
                    -90.0,
                    lat0_SO + 90,
                    lon0,
                    360.0,
                    cyclic_x=True,
                )
                spherical_cap.grid_metrics()
                print(
                    "spherical cap max/min latitude=",
                    spherical_cap.y.max(),
                    spherical_cap.y.min(),
                )
                print(
                    "spherical cap nj,ni=",
                    spherical_cap.y.shape[0] - 1,
                    spherical_cap.y.shape[1] - 1,
                )
                print("spherical cap starting longitude=", spherical_cap.x[0, 0])
                print("spherical cap ending longitude=", spherical_cap.x[0, -1])

                lamSC, phiSC = spherical_cap.x, spherical_cap.y
                dxSC, dySC, areaSC, angleSC = (
                    spherical_cap.dx,
                    spherical_cap.dy,
                    spherical_cap.area,
                    spherical_cap.angle_dx,
                )

                print(
                    "   CHECK_metrics_MIDAS: % errors in (area, lat arc, lon arc)",
                    metrics_error(dxSC, dySC, areaSC, Ni, phiSC[-1, 0], phiSC[0, 0]),
                )
        else:
            if not reproduce_MIDAS_grids:
                if match_dy:
                    # The SC grid should end at the same lattitude that SO starts.
                    # Then when we combine the two grids we should drop the x,y,dx,angle from one of the two.
                    # This way we get a continous dy and area.
                    lat0_SC = lat0_SO

                lamSC, phiSC = generate_displaced_pole_grid(
                    Ni, Nj_scap, lon0, lat0_SC, lon_dp, r_dp
                )
                angleSC = angle_x(lamSC, phiSC)
                if plotem:
                    displacedPoleCap_plot(
                        lamSC, phiSC, lon0, lon_dp, lat0_SO, stride=int(refineR * 10)
                    )
                    displacedPoleCap_plot2(
                        lamSC, phiSC, lon0, lon_dp, lat0_SO, stride=int(refineR * 10)
                    )

                # Quadrature metrics using great circle approximations for the h's
                dxSC = -np.ones([lamSC.shape[0], lamSC.shape[1] - 1])
                dySC = -np.ones([lamSC.shape[0] - 1, lamSC.shape[1]])
                areaSC = -np.ones([lamSC.shape[0] - 1, lamSC.shape[1] - 1])
                if calculate_metrics:
                    dxSC, dySC, areaSC = displacedPoleCap_metrics_quad(
                        4, Ni, Nj_scap, lon0, lat0_SC, lon_dp, r_dp
                    )
                    poles_i = int(Ni * np.mod(lon_dp - lon0, 360) / 360.0)
                    print(
                        "   CHECK_metrics_hquad: % errors in (area, lat arc, lon arc)",
                        metrics_error(
                            dxSC,
                            dySC,
                            areaSC,
                            Ni,
                            lat1=lat0_SC,
                            lat2=-90.0,
                            displaced_pole=poles_i,
                            excluded_fraction=doughnut,
                        ),
                    )
                # Compare with MIDAS metrics
                # dxSC_MIDAS,dySC_MIDAS,areaSC_MIDAS = generate_grid_metrics_MIDAS(lamSC,phiSC)
                # print("   CHECK_metrics_MIDAS: % errors in (area, lat arc, lon arc)", metrics_error(dxSC,dySC,areaSC,Ni,phiSC[-1,0],phiSC[0,0]))

                # Cut the unused portion of the grid
                # Choose the doughnut factor to keep the number of j's the same as in existing OM4p25 grid
                if doughnut != 0.0:
                    jmin = np.ceil(doughnut * Nj_scap)
                    jmin = jmin + np.mod(jmin, 2)
                    jmint = int(jmin)
                    lamSC = lamSC[jmint:, :]
                    phiSC = phiSC[jmint:, :]
                    dxSC = dxSC[jmint:, :]
                    dySC = dySC[jmint:, :]
                    areaSC = areaSC[jmint:, :]
                    angleSC = angleSC[jmint:, :]

                if phiSC.shape[0] % 2 == 0 and ensure_nj_even:
                    print(
                        "   The number of j's is not even. Fixing this by cutting one row at south."
                    )
                    lamSC = np.delete(lamSC, 0, 0)
                    phiSC = np.delete(phiSC, 0, 0)
                    dxSC = np.delete(dxSC, 0, 0)
                    dySC = np.delete(dySC, 0, 0)
                    areaSC = np.delete(areaSC, 0, 0)
                    angleSC = np.delete(angleSC, 0, 0)

                print("   number of js=", lamSC.shape[0])
            else:
                ny_scap = int(refineR * 40)
                r0_pole = 0.20
                lon0_pole = 100.0  # 180-longitude of the displaced pole
                doughnut = 0.12
                print(spherical.dy.shape)
                lenlat = 90.0 + spherical.y.min()
                dy0 = spherical.dy[0, 0] * r0_pole
                x = spherical.x[0, :]
                y = np.linspace(-90.0, 0.5 * (lat0_SO - 90.0), ny_scap / 8)
                y = np.concatenate(
                    (y, np.linspace(y.max(), lat0_SO, 7 * ny_scap / 8 + 1))
                )
                X, Y = np.meshgrid(x, y)
                antarctic_cap = supergrid(
                    xdat=X,
                    ydat=Y,
                    axis_units="degrees",
                    displace_pole=True,
                    r0_pole=r0_pole,
                    lon0_pole=lon0_pole,
                    doughnut=doughnut,
                )
                antarctic_cap.grid_metrics()
                print(
                    "generated a southern cap of size (ny,nx)= ",
                    antarctic_cap.y.shape[0] - 1,
                    antarctic_cap.y.shape[1] - 1,
                )

                lamSC, phiSC = antarctic_cap.x, antarctic_cap.y
                dxSC, dySC, areaSC, angleSC = (
                    antarctic_cap.dx,
                    antarctic_cap.dy,
                    antarctic_cap.area,
                    antarctic_cap.angle_dx,
                )

                print(
                    "   CHECK_metrics_MIDAS: % errors in (area, lat arc, lon arc)",
                    metrics_error(dxSC, dySC, areaSC, Ni, phiSC[-1, 0], phiSC[0, 0]),
                )

        if not "all" in grids:  # if only "sc" was requested cut it according to args
            # Cut the grid at south according to the options!
            # We may need to cut the whole SC grid and some of the SO
            cut = False
            jcut = 0
            cats = 0
            if south_cutoff_row > 0:
                cut = True
                jcut = south_cutoff_row - 1
            elif south_cutoff_ang > -90:
                cut = True
                jcut = 1 + np.nonzero(phiSC[:, 0] < south_cutoff_ang)[0][-1]

            if cut:
                print("   SC: shape[0], jcut", lamSC.shape[0], jcut)
                if jcut < lamSC.shape[0]:  # only SC needs to be cut
                    if (phiSC.shape[0] - jcut) % 2 == 0 and ensure_nj_even:
                        # if((areaSC.shape[0]-jcut-1)%2 == 0 and ensure_nj_even):
                        print(
                            "   SC: The number of j's is not even. Fixing this by cutting one row at south."
                        )
                        jcut = jcut + 1
                    print("   Cutting SC grid rows 0 to ", jcut)
                    lamSC = lamSC[jcut:, :]
                    phiSC = phiSC[jcut:, :]
                    dxSC = dxSC[jcut:, :]
                    dySC = dySC[jcut:, :]
                    areaSC = areaSC[jcut:, :]
                    angleSC = angleSC[jcut:, :]

        if write_subgrid_files:
            write_nc(
                lamSC,
                phiSC,
                dxSC,
                dySC,
                areaSC,
                angleSC,
                axis_units="degrees",
                fnam=gridfilename + "SC.nc",
                description=desc,
                history=hist,
                source=source,
                debug=debug,
            )

    if "all" in grids:
        # Concatenate to generate the whole grid
        # Start from displaced southern cap and join the southern ocean grid
        print("Stitching the grids together...")

        # Note that x,y,dx,angle_dx have a shape[0]=nyp1 and should be cut by one (total 3) for the merged variables
        # to have the right shape.
        # But y and area have a shape[0]=ny and should not be cut.
        # Niki: This cut can be done in a few ambigous ways:
        #    1.  MIDAS way (above)
        #    2.  this way  : x1=np.concatenate((lamSC[:-1,:],lamSO),axis=0)
        #                    y1=np.concatenate((phiSC[:-1,:],phiSO),axis=0)
        #                    dx1=np.concatenate((dxSC[:-1,:],dxSO),axis=0)
        #                    angle1=np.concatenate((angleSC[:-1,:],angleSO),axis=0)
        #                    #
        #                    dy1=np.concatenate((dySC,dySO),axis=0)
        #                    area1=np.concatenate((areaSC,areaSO),axis=0)
        #     3.  at the very end by restricting to x3[2:,:], ...
        #      Which way is "correct"?
        #      If the sub-grids are disjoint, 1 and 2 introduce a jump in y1 at the joint,
        #      as a result y1 and dy1 may become inconsistent?
        #
        hasSC = True
        if match_dy:
            # Cut the grid at south according to the options!
            # We may need to cut the whole SC grid and some of the SO
            cut = False
            jcut = 0
            cats = 0
            if south_cutoff_row > 0:
                cut = True
                jcut = south_cutoff_row - 1
            elif south_cutoff_ang > -90:
                cut = True
                jcut = 1 + np.nonzero(phiSC[:, 0] < south_cutoff_ang)[0][-1]

            if cut:
                print("   SC: shape[0], jcut", lamSC.shape[0], jcut)
                if jcut < lamSC.shape[0]:  # only SC needs to be cut
                    if (phiSC.shape[0] - jcut) % 2 == 0 and ensure_nj_even:
                        # if((areaSC.shape[0]-jcut-1)%2 == 0 and ensure_nj_even):
                        print(
                            "   SC: The number of j's is not even. Fixing this by cutting one row at south."
                        )
                        jcut = jcut + 1
                    print("   Cutting SC grid rows 0 to ", jcut)
                    lamSC = lamSC[jcut:, :]
                    phiSC = phiSC[jcut:, :]
                    dxSC = dxSC[jcut:, :]
                    dySC = dySC[jcut:, :]
                    areaSC = areaSC[jcut:, :]
                    angleSC = angleSC[jcut:, :]
                else:
                    print("   Whole SC and some of SO need to be cut!")
                    print("   SO: shape[0], jcut", lamSO.shape[0], jcut)
                    hasSC = False
                    jcut_SO = jcut - lamSC.shape[0]
                    #                    jcut_SO = max(jcut-lamSC.shape[0], 1 + np.nonzero(phiSO[:,0] < south_cutoff_ang)[0][-1])
                    if (areaSO.shape[0] - jcut_SO - 1) % 2 == 0 and ensure_nj_even:
                        print(
                            "   SO: The number of j's is not even. Fixing this by cutting one row at south."
                        )
                        jcut_SO = jcut_SO + 1
                    print("   No SC grid remained. Cutting SO grid rows 0 to ", jcut_SO)
                    lamSO = lamSO[jcut_SO:, :]
                    phiSO = phiSO[jcut_SO:, :]
                    dxSO = dxSO[jcut_SO:, :]
                    dySO = dySO[jcut_SO:, :]
                    areaSO = areaSO[jcut_SO:, :]
                    angleSO = angleSO[jcut_SO:, :]

            if hasSC:
                x1 = np.concatenate((lamSC[:-1, :], lamSO), axis=0)
                y1 = np.concatenate((phiSC[:-1, :], phiSO), axis=0)
                dx1 = np.concatenate((dxSC[:-1, :], dxSO), axis=0)
                angle1 = np.concatenate((angleSC[:-1, :], angleSO), axis=0)
                #
                dy1 = np.concatenate((dySC, dySO), axis=0)
                area1 = np.concatenate((areaSC, areaSO), axis=0)
                cats = cats + 1
            else:  # if the whole SC was cut
                x1 = lamSO
                y1 = phiSO
                dx1 = dxSO
                dy1 = dySO
                area1 = areaSO
                angle1 = angleSO

            # Join the Mercator grid
            x2 = np.concatenate((x1[:-1, :], lamMerc), axis=0)
            y2 = np.concatenate((y1[:-1, :], phiMerc), axis=0)
            dx2 = np.concatenate((dx1[:-1, :], dxMerc), axis=0)
            angle2 = np.concatenate((angle1[:-1, :], angleMerc), axis=0)
            #
            dy2 = np.concatenate((dy1, dyMerc), axis=0)
            area2 = np.concatenate((area1, areaMerc), axis=0)
            cats = cats + 1
            # Join the norhern bipolar cap grid
            x3 = np.concatenate((x2[:-1, :], lamBP), axis=0)
            y3 = np.concatenate((y2[:-1, :], phiBP), axis=0)
            dx3 = np.concatenate((dx2[:-1, :], dxBP), axis=0)
            angle3 = np.concatenate((angle2[:-1, :], angleBP), axis=0)
            #
            dy3 = np.concatenate((dy2, dyBP), axis=0)
            area3 = np.concatenate((area2, areaBP), axis=0)
            dy3_ = np.roll(y3[:, Ni // 4], shift=-1, axis=0) - y3[:, Ni // 4]
            if np.any(dy3_ == 0):
                raise Exception(
                    "lattitude array has repeated values along symmetry meridian!"
                )

        else:
            x1 = np.concatenate((lamSC, lamSO[1:, :]), axis=0)
            y1 = np.concatenate((phiSC, phiSO[1:, :]), axis=0)
            dx1 = np.concatenate((dxSC, dxSO[1:, :]), axis=0)
            dy1 = np.concatenate((dySC, dySO), axis=0)
            area1 = np.concatenate((areaSC, areaSO), axis=0)
            angle1 = np.concatenate((angleSC[:-1, :], angleSO[:-1, :]), axis=0)
            # Join the Mercator grid
            x2 = np.concatenate((x1, lamMerc[1:, :]), axis=0)
            y2 = np.concatenate((y1, phiMerc[1:, :]), axis=0)
            dx2 = np.concatenate((dx1, dxMerc[1:, :]), axis=0)
            dy2 = np.concatenate((dy1, dyMerc), axis=0)
            area2 = np.concatenate((area1, areaMerc), axis=0)
            angle2 = np.concatenate((angle1, angleMerc[:-1, :]), axis=0)
            # Join the norhern bipolar cap grid
            x3 = np.concatenate((x2, lamBP[1:, :]), axis=0)
            y3 = np.concatenate((y2, phiBP[1:, :]), axis=0)
            dx3 = np.concatenate((dx2, dxBP[1:, :]), axis=0)
            dy3 = np.concatenate((dy2, dyBP), axis=0)
            area3 = np.concatenate((area2, areaBP), axis=0)
            angle3 = np.concatenate((angle2, angleBP), axis=0)

            if south_cutoff_row > 0:
                jcut = south_cutoff_row - 1
                print("cutting grid rows 0 to ", jcut)
                x3 = x3[jcut:, :]
                y3 = y3[jcut:, :]
                dx3 = dx3[jcut:, :]
                dy3 = dy3[jcut:, :]
                area3 = area3[jcut:, :]
                angle3 = angle3[jcut:, :]

            if south_cutoff_ang > -90:
                jcut = 1 + np.nonzero(y3[:, 0] < south_cutoff_ang)[0][-1]
                print("cutting grid below ", south_cutoff_ang, jcut)
                x3 = x3[jcut:, :]
                y3 = y3[jcut:, :]
                dx3 = dx3[jcut:, :]
                dy3 = dy3[jcut:, :]
                area3 = area3[jcut:, :]
                angle3 = angle3[jcut:, :]

        dy3_ = np.roll(y3[:, Ni // 4], shift=-1, axis=0) - y3[:, Ni // 4]
        if np.any(dy3_ == 0):
            print(
                "WARNING: lattitude array has repeated values along symmetry meridian! Try option --match_dy"
            )

        if write_subgrid_files:
            if hasSC:
                write_nc(
                    lamSC,
                    phiSC,
                    dxSC,
                    dySC,
                    areaSC,
                    angleSC,
                    axis_units="degrees",
                    fnam=gridfilename + "SC.nc",
                    description=desc,
                    history=hist,
                    source=source,
                    debug=debug,
                )
            else:
                print(
                    "There remained no South Pole cap grid because of the number of rows cut= ",
                    jcut,
                    lamSC.shape[0],
                )
        if plotem:
            displacedPoleCap_plot(
                lamSC, phiSC, lon0, lon_dp, lat0_SO, stride=int(refineR * 10)
            )
            displacedPoleCap_plot2(
                lamSC,
                phiSC,
                lon0,
                lon_dp,
                lat0_SO,
                stride=int(refineR * 10),
                block=True,
            )

        # write the whole grid file
        desc = (
            desc
            + "It consists of a Mercator grid spanning "
            + str(phiMerc[0, 0])
            + " to "
            + str(phiMerc[-1, 0])
            + " degrees, flanked by a bipolar northern cap and a regular lat-lon grid spanning "
            + str(phiMerc[0, 0])
            + " to "
            + str(lat0_SO)
            + " degrees. "
        )

        desc = desc + "It is capped by a "
        if r_dp != 0.0:
            desc = desc + "displaced pole "
        else:
            desc = desc + "regular "
        desc = desc + "grid south of " + str(lat0_SO) + " degrees."

        if south_cutoff_ang > -90:
            desc = desc + " It is cut south of " + str(south_cutoff_ang) + " degrees."

        if south_cutoff_row > 0:
            desc = (
                desc
                + " The first "
                + str(south_cutoff_row)
                + " rows at south are deleted."
            )

        # Ensure that the equator (y=0) is still a u-point
        equator = 0.0
        equator_index = np.searchsorted(y3[:, Ni // 4], equator)
        if equator_index == 0:
            raise Exception("   Ooops: Equator is not in the grid")
        else:
            print("   Equator is at j=", equator_index)
        # Ensure that the equator (y=0) is a u-point
        if equator_index % 2 == 0:
            raise Exception(
                "Ooops: Equator is not going to be a u-point. Use option --south_cutoff_row to one more or on less row from south."
            )
        if y3.shape[0] % 2 == 0:
            raise Exception(
                "Ooops: The number of j's in the supergrid is not even. Use option --south_cutoff_row to one more or on less row from south."
            )
        if not reproduce_MIDAS_grids:
            write_nc(
                x3,
                y3,
                dx3,
                dy3,
                area3,
                angle3,
                axis_units="degrees",
                fnam=gridfilename,
                description=desc,
                history=hist,
                source=source,
                no_changing_meta=no_changing_meta,
                debug=debug,
            )
        else:
            write_nc(
                x3,
                y3,
                dx3,
                dy3,
                area3,
                angle3,
                axis_units="degrees",
                fnam=gridfilename,
                format="NETCDF3_CLASSIC",
                description=desc,
                history=hist,
                source=source,
                no_changing_meta=no_changing_meta,
                debug=debug,
            )

        print("Wrote the whole grid to file ", gridfilename)

        # Visualization
        if plotem:
            plot_mesh_in_xyz(
                x2, y2, stride=30, upperlat=-40, title="Grid south of -40 degrees"
            )
            plot_mesh_in_xyz(
                x3, y3, stride=30, lowerlat=40, title="Grid north of 40 degrees"
            )


if __name__ == "__main__":
    main(sys.argv[1:])
