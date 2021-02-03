"""
Module containing code to work with ocean transports at SAMBA array

"""

import numpy as np
import copy

from metric import constants
from metric import utils
from metric.geostrophy \
   import calc_dh, calc_vgeo, calc_vgeo_td, calc_vgeo_endpoint, \
          calc_vgeo_td_endpoint, update_georef, merge_vgeo_and_v, \
          calc_ek
from metric.samba import output


class SAMBA_Transports(object):
    """ Class to interface with volume and heat transport diagnostics """

    def __init__(self, v, t_on_v, s_on_v, ssh_on_v, dynht, rho, minind, maxind, sref=35.17):
        """ Initialize with velocity and temperature sections """

        # Initialize data
        self.name = v.name
        self.v = v.data[:,:,minind:maxind]
        self.t = t_on_v.data[:,:,minind:maxind]
        self.s = s_on_v.data[:,:,minind:maxind]
        if ssh_on_v is not None:
          self.ssh = ssh_on_v.data[:,minind:maxind]
        self.dh = dynht[:,:,minind:maxind+1]
        self.rho = rho[:,:,minind:maxind+1]
        self.rhocp = constants.RHO_REF * constants.CP
        self.sref = sref
        self.x = v.x[minind:maxind]
        self.xbounds = v.xbounds[minind:maxind+1]
        self.y = v.y[minind:maxind]
        self.z = v.z
        self.zbounds = v.zbounds
        self.dates = v.dates
        self.dz = v.dz
        self.dz_as_data = v.dz_as_data[:,:,minind:maxind]
        self.dx = v.cell_widths[minind:maxind]
        self.dx_as_data = v.cell_widths_as_data[:,:,minind:maxind]
        self.da = self.dx_as_data * self.dz_as_data

        # Set null values for property attributes
        self._avg_t = None
        self._avg_s = None
        self._avg_v = None
        self._v_no_net = None
        self._net_transport = None
        self._zonal_avg_v = None
        self._zonal_avg_v_no_net = None
        self._zonal_anom_v = None
        self._zonal_sum_v_no_net = None
        self._zonal_sum_v = None
        self._zonal_avg_t = None
        self._zonal_anom_t = None
        self._zonal_avg_s = None
        self._zonal_anom_s = None
        self._streamfunction = None
        self._oht_total = None
        self._oht_by_net = None
        self._oht_by_horizontal = None
        self._oht_by_overturning = None
        self._oft_total = None
        self._oft_by_net = None
        self._oft_by_horizontal = None
        self._oft_by_overturning = None

    def section_avg(self, data, total=False):
        """ Return avg across whole section """
        if total:
            return (data * self.da).sum(axis=(1,2))
        else:
            return (data * self.da).sum(axis=(1,2)) / self.da.sum(axis=(1,2))

    def zonal_avg(self, data, total=False):
        """ Return zonal avg across section """
        if total:
            return (data * self.dx_as_data).sum(axis=2)
        else:
            return (data * self.dx_as_data).sum(axis=2) / self.dx_as_data.sum(axis=2)

    @property
    def streamfunction(self):
        """ Return contribution to basin-wide overturning streamfunction """
        if self._streamfunction is None:
            self._streamfunction = np.cumsum(
                (self.v.filled(0) * self.da.filled(0)).sum(axis=2), axis=1)
        return self._streamfunction

    @property
    def avg_t(self):
        """ Return section average temperature """
        if self._avg_t is None:
            self._avg_t = self.section_avg(self.t)
        return self._avg_t

    @property
    def avg_s(self):
        """ Return section average salinity """
        if self._avg_s is None:
            self._avg_s = self.section_avg(self.s)
        return self._avg_s

    @property
    def avg_v(self):
        """ Return section average velocity """
        if self._avg_v is None:
            self._avg_v = self.section_avg(self.v)
        return self._avg_v

    @property
    def v_no_net(self):
        """ Return velocity after removing net transport through section """
        if self._v_no_net is None:
            self._v_no_net = self.v - self.avg_v[:,np.newaxis,np.newaxis]
        return self._v_no_net

    @property
    def net_transport(self):
        """ Return net transport through section """
        if self._net_transport is None:
            self._net_transport = self.section_avg(self.v, total=True)
        return self._net_transport

    @property
    def zonal_avg_v_no_net(self):
        """ Return zonal mean of v_no_net """
        if self._zonal_avg_v_no_net is None:
            self._zonal_avg_v_no_net = self.zonal_avg(self.v_no_net)
        return self._zonal_avg_v_no_net

    @property
    def zonal_avg_v(self):
        """ Return zonal mean of v """
        if self._zonal_avg_v is None:
            self._zonal_avg_v = self.zonal_avg(self.v)
        return self._zonal_avg_v

    @property
    def zonal_avg_t(self):
        """ Return zonal mean temperature profile """
        if self._zonal_avg_t is None:
            self._zonal_avg_t = self.zonal_avg(self.t)
        return self._zonal_avg_t

    @property
    def zonal_avg_s(self):
        """ Return zonal mean salinity profile """
        if self._zonal_avg_s is None:
            self._zonal_avg_s = self.zonal_avg(self.s)
        return self._zonal_avg_s

    @property
    def zonal_anom_v(self):
        """ Return velocity anomalies relative to zonal mean profile """
        if self._zonal_anom_v is None:
            self._zonal_anom_v = self.v_no_net - self.zonal_avg_v_no_net[:,:,np.newaxis]
        return self._zonal_anom_v

    @property
    def zonal_anom_t(self):
        """ Return temperature anomalies relative to zonal mean profile """
        if self._zonal_anom_t is None:
            self._zonal_anom_t = self.t - self.zonal_avg_t[:,:,np.newaxis]
        return self._zonal_anom_t

    @property
    def zonal_anom_s(self):
        """ Return salinity anomalies relative to zonal mean profile """
        if self._zonal_anom_s is None:
            self._zonal_anom_s = self.s - self.zonal_avg_s[:,:,np.newaxis]
        return self._zonal_anom_s

    @property
    def zonal_sum_v(self):
        """ Return zonal integral of v """
        if self._zonal_sum_v is None:
            self._zonal_sum_v = self.zonal_avg(self.v, total=True)
        return self._zonal_sum_v

    @property
    def zonal_sum_v_no_net(self):
        """ Return zonal integral of v_no_net """
        if self._zonal_sum_v_no_net is None:
            self._zonal_sum_v_no_net = self.zonal_avg(self.v_no_net, total=True)
        return self._zonal_sum_v_no_net

    @property
    def oht_by_net(self):
        """ Return heat transport by net transport through section """
        if self._oht_by_net is None:
            self._oht_by_net = self.net_transport * self.avg_t * self.rhocp
        return self._oht_by_net

    @property
    def oht_total(self):
        """ Return total heat transport through section """
        if self._oht_total is None:
            self._oht_total = self.section_avg(self.v * self.t, total=True) * self.rhocp
        return self._oht_total

    @property
    def oht_by_horizontal(self):
        """ Return heat transport by horizontal circulation """
        if self._oht_by_horizontal is None:
            self._oht_by_horizontal = self.section_avg(self.zonal_anom_v * self.zonal_anom_t,
                                                    total=True) * self.rhocp
        return self._oht_by_horizontal

    @property
    def oht_by_overturning(self):
        """ Return heat transport by local overturning circulation """
        if self._oht_by_overturning is None:
            self._oht_by_overturning = (self.zonal_sum_v_no_net * self.zonal_avg_t *
                                        self.dz[np.newaxis,:]).sum(axis=1) * self.rhocp
        return self._oht_by_overturning

    @property
    def oft_by_net(self):
        """ Return freshwater transport by net transport through section """
        if self._oft_by_net is None:
            self._oft_by_net = self.net_transport * (self.avg_s - self.sref) * (-1.0/self.sref)
        return self._oft_by_net

    @property
    def oft_total(self):
        """ Return total freshwater transport through section """
        if self._oft_total is None:
            self._oft_total =  self.section_avg(self.v * (self.s - self.sref), total=True) * (-1.0/self.sref)
        return self._oft_total

    @property
    def oft_by_horizontal(self):
        """ Return freshwater transport by horizontal circulation """
        if self._oft_by_horizontal is None:
            self._oft_by_horizontal = self.section_avg(self.zonal_anom_v * self.zonal_anom_s,
                                                    total=True) * (-1.0/self.sref)
        return self._oft_by_horizontal

    @property
    def oft_by_overturning(self):
        """ Return freshwater transport by local overturning circulation """
        if self._oft_by_overturning is None:
            self._oft_by_overturning = (self.zonal_sum_v_no_net * self.zonal_avg_s *
                                        self.dz[np.newaxis,:]).sum(axis=1) * (-1.0/self.sref)
        return self._oft_by_overturning


def calc_transports_from_sections(config, v, tau, t_on_v, s_on_v, ssh_on_v):
    """
    High-level routine to call transport calculations and return
    integrated transports on SAMBA section as netcdf object

    """
    # Extract sub-section boundaries
    wbw_minlon = config.getfloat('options','wbw_minlon') # Minimum Longitude of WBW
    int_minlon = config.getfloat('options','int_minlon') # Longitude of WBW/interior boundary
    int_maxlon = config.getfloat('options','int_maxlon') # Longitude of interior/EBW boundary
    ebw_maxlon = config.getfloat('options','ebw_maxlon') # Maximum Longitude of EBW

    # Get salinity reference used for freshwater transports
    sref = config.getfloat('options','reference_salinity')

    # Get indices for sub-sections
    wbwmin, wbwmax = utils.get_indrange(v.x, wbw_minlon, int_minlon)  # WBW
    intmin, intmax = utils.get_indrange(v.x, int_minlon, int_maxlon) # Gyre interior
    ebwmin, ebwmax = utils.get_indrange(v.x, int_maxlon, ebw_maxlon) # EBW

    # Calculate dynamic heights
    eos = config.get('options','eos')
    dh, rho = calc_dh(t_on_v, s_on_v, eos, get_rho=True)

    # Calculate geostrophic transports
    if config.getboolean('options', 'td_geo'):
      if config.getboolean('options', 'endpoint'):
        vgeo = calc_vgeo_td_endpoint(v, rho, ssh_on_v, intmin, intmax-1)
      else:
        vgeo = calc_vgeo_td(v, rho, ssh_on_v)
    else:
      georef = config.getfloat('options', 'georef_level')
      if config.getboolean('options', 'endpoint'):
        vgeo = calc_vgeo_endpoint(v, dh, intmin, intmax-1, georef=georef)
      else:  
        vgeo = calc_vgeo(v, dh, georef=georef)
      # Optionally reference geostrophic transports to model velocities
      if config.has_option('options', 'vref_level'):
        vref_level = config.getfloat('options', 'vref_level')
        vgeo = update_georef(vgeo, v, vref_level)

    # Calculate Ekman velocities
    ek_level = config.getfloat('options','ekman_depth')
    if config.has_option('options','ek_profile_type'):
        ek_profile = config.get('options','ek_profile_type')
    else:
        ek_profile = 'uniform'

    ek = calc_ek(v, tau, int_minlon, int_maxlon, ek_level, profile=ek_profile)

    # Use model velocities in WBW and EBW regions
    vgeo = merge_vgeo_and_v(vgeo, v, wbw_minlon, int_minlon)
    vgeo = merge_vgeo_and_v(vgeo, v, int_maxlon, ebw_maxlon)

    # Add ekman to geostrophic transports for combined samba velocities
    vsamba = copy.deepcopy(vgeo)
    vsamba.data = vgeo.data + ek.data

    # Get volume and heat transports on each (sub-)section
    # Western-boundary wedge transports
    wbw_trans = SAMBA_Transports(vgeo, t_on_v, s_on_v, ssh_on_v, dh.bounds_data, rho.bounds_data, wbwmin, wbwmax, sref=sref)
    # Gyre interior transports
    int_trans = SAMBA_Transports(vgeo, t_on_v, s_on_v, ssh_on_v, dh.bounds_data, rho.bounds_data, intmin, intmax, sref=sref)
    # Gyre interior transports using model velocities
    int_trans_mod = SAMBA_Transports(v, t_on_v, s_on_v, ssh_on_v, dh.bounds_data, rho.bounds_data, intmin, intmax, sref=sref)
    # Eastern-boundary wedge transports
    ebw_trans = SAMBA_Transports(vgeo, t_on_v, s_on_v, ssh_on_v, dh.bounds_data, rho.bounds_data, ebwmin, ebwmax, sref=sref)
    # Ekman transports
    ek_trans = SAMBA_Transports(ek, t_on_v, s_on_v, ssh_on_v, dh.bounds_data, rho.bounds_data, intmin, intmax, sref=sref)
    # Total section transports using SAMBA approximation
    samba_trans = SAMBA_Transports(vsamba, t_on_v, s_on_v, ssh_on_v, dh.bounds_data, rho.bounds_data, wbwmin, ebwmax, sref=sref)
    # Total section transports using model velocities
    model_trans = SAMBA_Transports(v, t_on_v, s_on_v, ssh_on_v, dh.bounds_data, rho.bounds_data, wbwmin, ebwmax, sref=sref)

    # Create netcdf object for output/plotting
    trans = output.create_netcdf(config, samba_trans, model_trans, wbw_trans,
                                 ebw_trans, int_trans, int_trans_mod, ek_trans)

    return trans

