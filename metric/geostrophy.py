"""
Functions for calculating geostrophic currents.

"""

import numpy as np
import copy

from . import constants
from . import utils



def eos_insitu(t, s, p):
    """
    Returns in situ density of seawater as calculated by the NEMO
    routine eos_insitu.f90. Computes the density referenced to
    a specified depth/pressure from potential temperature and salinity
    using the Jackett and McDougall (1994) equation of state.

    """
    # Convert to double precision
    ptem = np.double(t)    # potential temperature (celcius)
    psal = np.double(s)    # salintiy (psu)
    depth = np.double(p)   # pressure (decibar) = depth (m)
    rau0 = np.double(1035) # volumic mass of reference (kg/m3)

    # Read into eos_insitu.f90 varnames
    zrau0r = 1.e0 / rau0
    zt = ptem
    zs = psal
    zh = depth
    zsr= np.sqrt(np.abs(psal))   # square root salinity

    # compute volumic mass pure water at atm pressure
    zr1 = ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt-9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594

    # seawater volumic mass atm pressure
    zr2 = ( ( ( 5.3875e-9*zt-8.2467e-7 ) *zt+7.6438e-5 ) *zt-4.0899e-3 ) *zt+0.824493
    zr3 = ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3
    zr4 = 4.8314e-4

    #  potential volumic mass (reference to the surface)
    zrhop = ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1

    # add the compression terms
    ze = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6
    zbw = (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4
    zb = zbw + ze * zs
    zd = -2.042967e-2
    zc =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896
    zaw = ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859 ) *zt - 4.721788
    za = ( zd*zsr + zc ) *zs + zaw
    zb1 =   (-0.1909078*zt+7.390729 ) *zt-55.87545
    za1 = ( ( 2.326469e-3*zt+1.553190)*zt-65.00517 ) *zt+1044.077
    zkw = ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638 ) *zt + 2098.925 ) *zt+190925.6
    zk0 = ( zb1*zsr + za1 )*zs + zkw

    # Caculate density
    prd = (  zrhop / (  1.0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - rau0  ) * zrau0r
    rho = (prd*rau0) + rau0

    return rho


def calc_dh(t_on_v, s_on_v, get_rho=True):
    """
    Return ZonalSections containing dynamic heights calculated from
    from temperature and salinity interpolated onto velocity boundaries.

    """
    # Calculate in situ density at bounds
    rho = copy.deepcopy(t_on_v)
    rho.data = None # Density not needed at v mid-points
    rho.bounds_data = eos_insitu(t_on_v.bounds_data, s_on_v.bounds_data,
                                 t_on_v.z_as_bounds_data)

    # Calculate dynamic height relative to a reference level
    dh = copy.deepcopy(rho)
    rho_anom = (rho.bounds_data - constants.RHO_REF) / constants.RHO_REF
    # Depth axis reversed for integral from sea-floor.
    dh.bounds_data = constants.G * np.cumsum((rho_anom * rho.dz_as_bounds_data)[:,::-1,:],
                                axis=1)[:,::-1,:]

    if get_rho:
      return dh,rho
    else:
      return dh


def calc_vgeo(v, dh, georef=4950.):
    """
    Return ZonalSections containing geostrophic velocities
    relative to specified reference level.

    """
    vgeo = copy.deepcopy(v) # Copy velocity data structure

    for nprof in range(len(vgeo.x)): # Loop through profiles
        if not v.mask[:,:,nprof].all():

            # Extract depth and dynamic height profiles at bounds
            z = dh.z
            z1 = dh.z_as_bounds_data[:,:,nprof]
            z2 = dh.z_as_bounds_data[:,:,nprof+1]
            dh1 = dh.bounds_data[:,:,nprof]
            dh2 = dh.bounds_data[:,:,nprof+1]

            # Coriolis parameter at profile location
            corf = 2 * constants.ROT * np.sin(np.pi * (vgeo.y[nprof]/180.) )

            # cell width along section
            dx = vgeo.cell_widths[nprof]

            # Clip reference depth using ocean floor.
            maxz = np.min([z1.max(),z2.max()])
            zref = min(georef, maxz)

            # Adjust dh to new reference level
            zind = utils.find_nearest(z,zref)
            dh1 -= dh1[:,zind]
            dh2 -= dh2[:,zind]

            # Calculate geostrophic velocity
            vgeo_profile = (-1. / corf) * ( (dh2 - dh1) / dx )
            vgeo.data[:,:,nprof] = vgeo_profile

    return vgeo


def calc_vgeo_endpoint(v, dh, idx1, idx2, georef=4950.):
    """
    Return ZonalSections containing geostrophic velocities
    relative to specified reference level using only the endpoints.

    """
    vgeo = copy.deepcopy(v) # Copy velocity data structure

    if not (v.mask[:,:,idx1].all() | v.mask[:,:,idx2].all()):

        # Extract depth and dynamic height profiles at bounds
        z = dh.z
        z1 = dh.z_as_bounds_data[:,:,idx1]
        z2 = dh.z_as_bounds_data[:,:,idx2+1]
        dh1 = dh.bounds_data[:,:,idx1]
        dh2 = dh.bounds_data[:,:,idx2+1]

        # Coriolis parameter at profile location
        corf = 2 * constants.ROT * np.sin(np.pi * (vgeo.y[idx1:idx2+1].mean()/180.) )

        # cell width along section
        dx = vgeo.cell_widths[idx1:idx2+1].sum()

        # Clip reference depth using ocean floor.
        maxz = np.min([z1.max(),z2.max()])
        zref = min(georef, maxz)

        # Adjust dh to new reference level
        zind = utils.find_nearest(z,zref)
        dh1 -= dh1[:,zind]
        dh2 -= dh2[:,zind]

        # Calculate geostrophic velocity
        vgeo_profile = (-1. / corf) * ( (dh2 - dh1) / dx )
        vgeo.data[:,:,idx1:idx2+1] = vgeo_profile[:,:,None]

    return vgeo

def calc_vgeo_td(v, rho, ssh_on_v):
    """
    Return ZonalSections containing geostrophic velocities
    computed from surface to sea-floor.

    """
    vgeo = copy.deepcopy(v) # Copy velocity data structure

    press = copy.deepcopy(rho)
    press.bounds_data = np.cumsum((constants.G * rho.bounds_data * rho.dz_as_bounds_data), axis=1)

    for nprof in range(len(vgeo.x)): # Loop through profiles
        if not v.mask[:,:,nprof].all():

            # Extract pressure at bounds
            press1 = press.bounds_data[:,:,nprof]
            press2 = press.bounds_data[:,:,nprof+1]

            # Extract ssh at bounds
            ssh1 = ssh_on_v.bounds_data[:,nprof]
            ssh2 = ssh_on_v.bounds_data[:,nprof+1]

            # Extract density
            dens = 0.5 * (rho.bounds_data[:,:,nprof] + rho.bounds_data[:,:,nprof+1])

            # Coriolis parameter at profile location
            corf = 2 * constants.ROT * np.sin(np.pi * (vgeo.y[nprof]/180.) )

            # cell width along section
            dx = vgeo.cell_widths[nprof]

            # Calculate geostrophic velocity
            vgeo_barotropic = ((constants.G / corf) * ( (ssh2 - ssh1) / dx))
            vgeo_baroclinic = ((1 / (corf * dens)) * ( (press2 - press1) / dx))
            vgeo.data[:,:,nprof] = vgeo_barotropic[:,None] + vgeo_baroclinic

    return vgeo


def calc_vgeo_td_endpoint(v, rho, ssh_on_v, idx1, idx2):
    """
    Return ZonalSections containing geostrophic velocities
    computed from surface to sea-floor using only the endpoints.

    """
    vgeo = copy.deepcopy(v) # Copy velocity data structure

    press = copy.deepcopy(rho)
    press.bounds_data = np.cumsum((constants.G * rho.bounds_data * rho.dz_as_bounds_data), axis=1)

    if not (v.mask[:,:,idx1].all() | v.mask[:,:,idx2].all()):

            # Extract pressure at bounds
            press1 = press.bounds_data[:,:,idx1]
            press2 = press.bounds_data[:,:,idx2+1]

            # Extract ssh at bounds
            ssh1 = ssh_on_v.bounds_data[:,idx1]
            ssh2 = ssh_on_v.bounds_data[:,idx2+1]

            # Extract density
            dens = 0.5 * (rho.bounds_data[:,:,idx1] + rho.bounds_data[:,:,idx2+1])

            # Coriolis parameter at profile location
            corf = 2 * constants.ROT * np.sin(np.pi * (vgeo.y[idx1:idx2+1].mean()/180.) )

            # cell width along section
            dx = vgeo.cell_widths[idx1:idx2+1].sum()

            # Calculate geostrophic velocity
            vgeo_barotropic = ((constants.G / corf) * ( (ssh2 - ssh1) / dx))
            vgeo_baroclinic = ((1 / (corf * dens)) * ( (press2 - press1) / dx))
            vgeo.data[:,:,idx1:idx2+1] = vgeo_barotropic[:,None,None] + vgeo_baroclinic[:,:,None]

    return vgeo


def update_georef(vgeo, v, vref_level):
    """
    Return vgeo after updating geostrophic reference depth by constraining
    velocities in vgeo to match those in v at the specified depth.

    """
    vgeodat = vgeo.data.filled(0)
    vdat = v.data.filled(0)
    zind = utils.find_nearest(v.z,vref_level)
    vadj = np.ones_like(vgeodat) * (vdat[:,zind,:] - vgeodat[:,zind,:])
    vgeo.data = np.ma.MaskedArray(vgeo.data + vadj, mask=vgeo.mask)

    return vgeo


def calc_ek(v, tau, minlon, maxlon, ek_level, profile='uniform'):
    """ Return ZonalSections containing Ekman velocities """

    # Copy velocity data structure
    ek = copy.deepcopy(v)
    ek.data.data[:] = 0.0

    # Get indices for gyre interior
    intmin, intmax = utils.get_indrange(tau.x, minlon, maxlon)

    # Calculate depth-integrated Ekman transports
    dx = tau.cell_widths_as_data[:,intmin:intmax]
    lats = tau.y[intmin:intmax]
    taux = tau.data[:,intmin:intmax]
    corf = 2 * constants.ROT * np.sin(np.pi * (lats / 180.) )
    ek_trans = ((-1. *  taux / (corf * constants.RHO_REF)) * dx ).sum(axis=1)

    # Calculate velocities over ekman layer
    ek_minind, ek_maxind = utils.get_indrange(v.z, 0, ek_level)
    dz = v.dz_as_data[0,ek_minind:ek_maxind,intmin:intmax]
    dx = v.cell_widths_as_data[0,ek_minind:ek_maxind,intmin:intmax]

    if profile == 'uniform':
        # Use uniform Ekman velocities
        ek_area = (dx * dz).sum()
        ek.data[:,ek_minind:ek_maxind,intmin:intmax] = ek_trans[:,np.newaxis, np.newaxis] / ek_area
        ek.data = np.ma.MaskedArray(ek.data, mask=v.mask)
    elif profile == 'linear':
        # Use Ekman transport profile that linearly reduces to zero at z=zek
        zprof = ek.z[ek_minind:ek_maxind]
        dzprof = ek.dz[ek_minind:ek_maxind]
        zmax = dzprof.sum()
        vek = get_linear_profiles(ek_trans, zprof, dzprof, zmax) / dx[np.newaxis].sum(axis=2)
        ek.data[:,ek_minind:ek_maxind,intmin:intmax] = vek[:,:,np.newaxis]
        ek.data = np.ma.MaskedArray(ek.data, mask=v.mask)
    else:
        raise ValueError('Unrecognized ekman profile type')

    return ek


def get_linear_profiles(u_int, z, dz, zmax):
    """
    Return transport profile U_z that decreases linearly from z=0 to
    z=zmax and is constrained by u_int.

    u(z) = umax when z=0
    u(z) = 0 when z=zmax

    \int_{z=zmax}^{z=0} u(z) dz = u_int

    """

    u_max = 2 * u_int / zmax**2
    u_z = u_max[:,np.newaxis] * (zmax - z)[np.newaxis,:]
    scale = u_int / (u_z * dz[np.newaxis]).sum(axis=1) # Ensure conservation of integral

    return u_z * scale[:,np.newaxis]


def merge_vgeo_and_v(vgeo, v, minlon, maxlon):
    """ Return vgeo with velocities from v west of lonbnd """
    minind, maxind = utils.get_indrange(vgeo.x, minlon, maxlon)
    vgeo.data[:,:,minind:maxind] = v.data[:,:,minind:maxind]

    return vgeo


def section_integral(v, xmin, xmax):
    """ Section integral between x values """
    minind, maxind = utils.get_indrange(v.x, xmin, xmax)
    if v.surface_field:
        dx = v.cell_widths_as_data[:,minind:maxind]
        return np.sum(v.data[:,minind:maxind] * dx, axis=1)
    else:
        da = v.cell_widths_as_data[:,:,minind:maxind] * v.dz_as_data[:,:,minind:maxind]
        return np.sum(v.data[:,:,minind:maxind] * da , axis=(1,2))


def rapid_mass_balance(vgeo, ek, minlon, midlon, maxlon):
    """
    Return vgeo after applying RAPID-style mass-balance constraint
    as a barotropic velocity over geostrophic interior
    """

    # Calculate net transports
    fcwbw_tot = section_integral(vgeo, minlon, midlon)
    ek_tot = section_integral(ek, midlon, maxlon)
    int_tot = section_integral(vgeo, midlon, maxlon)
    net = int_tot + ek_tot + fcwbw_tot

    # Get cell dimensions in gyre interior
    minind,maxind = utils.get_indrange(vgeo.x, midlon, maxlon)
    dz = vgeo.dz_as_data[0]
    dx = vgeo.cell_widths_as_data[0]
    da = (dx[:,minind:maxind] * dz[:,minind:maxind])

    # Correct geostrophic transports in gyre interior
    corr = net / da.sum()
    vgeo.data[:,:,minind:maxind] = (vgeo.data[:,:,minind:maxind] -
                                    corr[:,np.newaxis,np.newaxis])

    return vgeo, net


def total_mass_balance(v):
    """ Apply mass balance evenly across entire section """
    da =  v.cell_widths_as_data * v.dz_as_data
    v.data = (v.data - ((v.data * da).sum(axis=(1,2)) /
                        da.sum(axis=(1,2)))[:,np.newaxis,np.newaxis])

    return v
