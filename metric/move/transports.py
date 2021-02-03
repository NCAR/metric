"""
Module containing code to work with ocean transports

"""

import numpy as np
import copy

from metric import constants
from metric import utils
from metric.geostrophy \
   import calc_dh, calc_vgeo, calc_vgeo_td, calc_vgeo_endpoint, \
          calc_vgeo_td_endpoint, update_georef, merge_vgeo_and_v
from metric.move import output


class MOVE_Transports(object):
    """ Class to interface with volume and heat transport diagnostics """

    def __init__(self, v, t_on_v, s_on_v, ssh_on_v, dynht, rho, xminind, xmaxind, zminind, zmaxind, sref=35.17):
        """ Initialize with velocity and temperature sections """

        # Initialize data
        self.name = v.name
        self.v = v.data[:,zminind:zmaxind,xminind:xmaxind]
        self.t = t_on_v.data[:,zminind:zmaxind,xminind:xmaxind]
        self.s = s_on_v.data[:,zminind:zmaxind,xminind:xmaxind]
        if ssh_on_v is not None:
          self.ssh = ssh_on_v.data[:,xminind:xmaxind]
        self.dh = dynht[:,zminind:zmaxind,xminind:xmaxind+1]
        self.rho = rho[:,zminind:zmaxind,xminind:xmaxind+1]
        self.rhocp = constants.RHO_REF * constants.CP
        self.sref = sref
        self.x = v.x[xminind:xmaxind]
        self.xbounds = v.xbounds[xminind:xmaxind+1]
        self.y = v.y[xminind:xmaxind]
        self.z = v.z[zminind:zmaxind]
        self.zbounds = v.zbounds[zminind:zmaxind+1]
        self.dates = v.dates
        self.dz = v.dz[zminind:zmaxind]
        self.dz_as_data = v.dz_as_data[:,zminind:zmaxind,xminind:xmaxind]
        self.dx = v.cell_widths[xminind:xmaxind]
        self.dx_as_data = v.cell_widths_as_data[:,zminind:zmaxind,xminind:xmaxind]
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


def calc_transports_from_sections(config, v, t_on_v, s_on_v, ssh_on_v):
    """
    High-level routine to call transport calculations and return
    integrated transports on MOVE section as netcdf object

    """

    # Extract upper/lower levels limit for the deep southward-flowing transport
    upper_level = config.getfloat('options','upper_level')   #  upper level of the NADW layer
    lower_level = config.getfloat('options','lower_level') # lower level coincides with the interface depth between AABW and NADW

    # Extract sub-section boundaries
    bdry_minlon = config.getfloat('options','bdry_minlon')   # Minimum longitude for boundary component
    bdry_maxlon = config.getfloat('options','bdry_maxlon')   # Longitude of boundary/internal component boundary
    int_maxlon = config.getfloat('options','int_maxlon')   # Maximum longitude for internal component boundary

    # Get salinity reference used for freshwater transports
    sref = config.getfloat('options','reference_salinity')

    # Get indices for sub-sections
    bdrymin, bdrymax = utils.get_indrange(v.x, bdry_minlon, bdry_maxlon)     # boundary component
    intmin,intmax = utils.get_indrange(v.x, bdry_maxlon, int_maxlon)         # boundary component

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
      georef = config.getfloat('options', 'georef_level') # deep level of no motion
      if config.getboolean('options', 'endpoint'):
        vgeo = calc_vgeo_endpoint(v, dh, intmin, intmax-1, georef=georef)
      else:
        vgeo = calc_vgeo(v, dh, georef=georef)
      # Optionally reference geostrophic transports to model velocities
      if config.has_option('options', 'vref_level'):
        vref_level = config.getfloat('options', 'vref_level')
        vgeo = update_georef(vgeo, v, vref_level)

    # Use model velocities for boundary component
    vgeo = merge_vgeo_and_v(vgeo, v, bdry_minlon, bdry_maxlon)

    # Get volume and heat transports on each (sub-)section
    zindmin,zindmax = utils.get_indrange(v.z,upper_level,lower_level)
    # Total section transports using MOVE approximation
    move_trans = MOVE_Transports(vgeo, t_on_v, s_on_v, ssh_on_v, dh.bounds_data, rho.bounds_data, bdrymin, intmax, zindmin, zindmax, sref=sref)
    # Total section transports using model velocities
    model_trans = MOVE_Transports(v, t_on_v, s_on_v, ssh_on_v, dh.bounds_data, rho.bounds_data, bdrymin, intmax, zindmin, zindmax, sref=sref)
    # boundary transports
    bdry_trans = MOVE_Transports(vgeo, t_on_v, s_on_v, ssh_on_v, dh.bounds_data, rho.bounds_data, bdrymin, bdrymax, zindmin, zindmax, sref=sref)
    # internal transports
    int_trans = MOVE_Transports(vgeo, t_on_v, s_on_v, ssh_on_v, dh.bounds_data, rho.bounds_data, intmin, intmax, zindmin, zindmax, sref=sref)
    # internal transports using model velocities
    int_model_trans = MOVE_Transports(v, t_on_v, s_on_v, ssh_on_v, dh.bounds_data, rho.bounds_data, intmin, intmax, zindmin, zindmax, sref=sref)

    # Create netcdf object for output/plotting
    trans = output.create_netcdf(config, move_trans,  model_trans, bdry_trans, int_trans, int_model_trans)

    return trans

