"""
Module containing code to generate netcdf object for output.

"""


from netCDF4 import Dataset, date2num, default_fillvals
import numpy as np
import os

from metric import utils


def open_ncfile(config, dates):
    """ Return handle for netcdf file """
    outdir = config.get('output', 'outdir')
    name = config.get('output', 'name')
    date_format = config.get('output', 'date_format')

    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    if config.getboolean('options', 'td_geo'):
      if config.getboolean('options', 'endpoint'):
        suffix='_natl_meridional_transports_at_34_5S_td_endpoint.nc'
      else:
        suffix='_natl_meridional_transports_at_34_5S_td.nc'
    else:
      if config.getboolean('options', 'endpoint'):
        suffix='_natl_meridional_transports_at_34_5S_endpoint.nc'
      else:
        suffix='_natl_meridional_transports_at_34_5S.nc'

    savef = utils.get_savename(
        outdir, name, dates, date_format,
        suffix=suffix)
    dataset = Dataset(savef, 'w', format='NETCDF4_CLASSIC')

    return dataset

def create_netcdf(config, samba_trans, model_trans, wbw_trans,
                  ebw_trans, int_trans, int_trans_mod, ek_trans):
    """
    Return diagnostics for comparison with SAMBA observations
    in netcdf object for plotting/output.

    """

    # Configuration options
    zind = samba_trans.streamfunction.mean(axis=0).argmax()
    wbw_minlon = config.getfloat('options','wbw_minlon')
    int_minlon = config.getfloat('options','int_minlon')
    int_maxlon = config.getfloat('options','int_maxlon')
    ebw_maxlon = config.getfloat('options','ebw_maxlon')
    georef = config.getfloat('options','georef_level')
    ek_level = config.getfloat('options','ekman_depth')
    eos = config.get('options','eos')

    try:
        vref_level =  config.getfloat('options','vref_level')
    except:
        vref_level = 'None'

    # Create netcdf file and add dimensions
    dataset = open_ncfile(config, samba_trans.dates)
    zdim = dataset.createDimension('z', samba_trans.z.size)
    zbdim = dataset.createDimension('zbounds', samba_trans.zbounds.size)
    xdim = dataset.createDimension('x', samba_trans.x.size)
    xbdim = dataset.createDimension('xbounds', samba_trans.xbounds.size)
    tdim = dataset.createDimension('time', None)

    # Create time coordinate
    time = dataset.createVariable('time',np.float64,(tdim.name,))
    time.units = 'hours since 0001-01-01 00:00:00.0'
    time.calendar = 'gregorian'
    time[:] = date2num(samba_trans.dates, time.units, calendar=time.calendar)

    # Create depth coordinate
    z = dataset.createVariable('z',np.float64,(zdim.name,))
    z.units = 'm'
    z[:] = samba_trans.z

    # Create depth coordinate
    zbounds = dataset.createVariable('zbounds',np.float64,(zbdim.name,))
    zbounds.units = 'm'
    zbounds[:] = samba_trans.zbounds

    # Create depth coordinate
    dz = dataset.createVariable('dz',np.float64,(zdim.name,))
    dz.units = 'm'
    dz[:] = samba_trans.dz

    # Create lon coordinate
    x = dataset.createVariable('x',np.float64,(xdim.name,))
    x.units = 'degrees_east'
    x[:] = samba_trans.x

    # Create lon coordinate
    dx = dataset.createVariable('dx',np.float64,(xdim.name,))
    dx.units = 'm'
    dx[:] = samba_trans.dx

    # Create lon_bounds coordinate
    xbounds = dataset.createVariable('xbounds',np.float64,(xbdim.name,))
    xbounds.units = 'degrees_east'
    xbounds[:] = samba_trans.xbounds

    # Global attributes
    if config.getboolean('options', 'td_geo'):
      dataset.geostrophic_method = 'top-down'
    else:
      dataset.geostrophic_method = 'bottom-up'
      dataset.geostrophic_reference_level = georef
      dataset.reference_to_model_velocity = vref_level
    if config.getboolean('options', 'endpoint'):
      dataset.geostrophic_computation = 'endpoint'
    dataset.eos = eos
    dataset.rhocp = samba_trans.rhocp
    dataset.ekman_level = ek_level
    dataset.contact = 'fredc.ucar.edu'
    dataset.code_reference = 'https://github.com/NCAR/metric'
    dataset.method_references = 'Kersale, M., and Coauthors, 2020: Highly variable upper and abyssal overturning cells in the South Atlantic, Science Advances, Vol. 6, no. 32, https://doi.org/10.1126/sciadv.aba7573'

    # SAMBA velocity cross section
    v = dataset.createVariable('v',np.float64,(tdim.name,zdim.name,xdim.name),fill_value=default_fillvals['f8'])
    v.units = 'm/s'
    v.comment = 'model velocity'
    v[:] = model_trans.v

    # SAMBA geo velocity cross section
    vgeo = dataset.createVariable('vgeo',np.float64,(tdim.name,zdim.name,xdim.name),fill_value=default_fillvals['f8'])
    vgeo.units = 'm/s'
    vgeo.comment = 'geostrophic velocity'
    vgeo[:] = samba_trans.v

    # SAMBA dh cross section
    dh = dataset.createVariable('dh',np.float64,(tdim.name,zdim.name,xbdim.name),fill_value=default_fillvals['f8'])
    dh.units = 'm2/s2'
    dh.comment = 'dynamic heights'
    dh[:] = samba_trans.dh

    # SAMBA rho cross section                                                                                                        
    rho = dataset.createVariable('rho',np.float64,(tdim.name,zdim.name,xbdim.name),fill_value=default_fillvals['f8'])               
    rho.units = 'kg/m3'                                                                                                             
    rho.comment = 'density'                                                                                                 
    rho[:] = samba_trans.rho

    # SAMBA temp cross section
    temp = dataset.createVariable('temp',np.float64,(tdim.name,zdim.name,xdim.name),fill_value=default_fillvals['f8'])
    temp.units = 'degC'
    temp.comment = 'model potential temperature'
    temp[:] = model_trans.t

    # SAMBA salinity cross section
    salt = dataset.createVariable('salt',np.float64,(tdim.name,zdim.name,xdim.name),fill_value=default_fillvals['f8'])
    salt.units = 'PSU'
    salt.comment = 'model salinity'
    salt[:] = model_trans.s

    # SAMBA SSH cross section
    if config.getboolean('options', 'td_geo'):
      ssh = dataset.createVariable('ssh',np.float64,(tdim.name,xdim.name),fill_value=default_fillvals['f8'])
      ssh.units = 'm'
      ssh.comment = 'model SSH'
      ssh[:] = model_trans.ssh

    # Basinwide potential temperature profile
    t_basin = dataset.createVariable('t_basin',np.float64,(tdim.name,zdim.name))
    t_basin.units = 'degC'
    t_basin.minimum_longitude = wbw_minlon
    t_basin.maximum_longitude = ebw_maxlon
    t_basin.comment = 'Basinwide zonal mean potential temperature profile'
    t_basin[:] = samba_trans.zonal_avg_t

    # Basinwide salinity profile
    s_basin = dataset.createVariable('s_basin',np.float64,(tdim.name,zdim.name))
    s_basin.units = 'PSU'
    s_basin.minimum_longitude = wbw_minlon
    s_basin.maximum_longitude = ebw_maxlon
    s_basin.comment = 'Basinwide zonal mean salinity profile'
    s_basin[:] = samba_trans.zonal_avg_s

    # Basinwide transport profile - SAMBA approx
    v_basin_samba = dataset.createVariable('v_basin_samba',np.float64,(tdim.name,zdim.name))
    v_basin_samba.units = 'Sv/m'
    v_basin_samba.minimum_longitude = wbw_minlon
    v_basin_samba.maximum_longitude = ebw_maxlon
    v_basin_samba.comment = 'Basinwide transport profile using SAMBA approximations'
    v_basin_samba[:] = samba_trans.zonal_sum_v / 1e6

    # Basinwide transport profile - model v
    v_basin_model = dataset.createVariable('v_basin_model',np.float64,(tdim.name,zdim.name))
    v_basin_model.units = 'Sv/m'
    v_basin_model.minimum_longitude = wbw_minlon
    v_basin_model.maximum_longitude = ebw_maxlon
    v_basin_model.comment = 'Basinwide transport profile using model velocities'
    v_basin_model[:] = model_trans.zonal_sum_v / 1e6

    # Ekman transport time series
    ekman = dataset.createVariable('ekman',np.float64,(tdim.name))
    ekman.units = 'Sv'
    ekman.minimum_longitude = int_minlon
    ekman.maximum_longitude = int_maxlon
    ekman.comment = 'Ekman transport time series (streamfunction at 1000m)'
    ekman[:] = ek_trans.streamfunction[:,zind] / 1e6

    # Gyre interior transport time series
    geoint = dataset.createVariable('geoint',np.float64,(tdim.name))
    geoint.units = 'Sv'
    geoint.minimum_longitude = int_minlon
    geoint.maximum_longitude = int_maxlon
    geoint.comment = 'Geostrophic interior transport time series (streamfunction at 1000m).'
    geoint[:] = int_trans.streamfunction[:,zind] / 1e6

   # Gyre interior transport time series
    int_mod = dataset.createVariable('int_mod',np.float64,(tdim.name))
    int_mod.units = 'Sv'
    int_mod.minimum_longitude = int_minlon
    int_mod.maximum_longitude = int_maxlon
    int_mod.comment = 'Model interior transport time series (streamfunction at 1000m).'
    int_mod[:] = int_trans_mod.streamfunction[:,zind] / 1e6

    # Western-boundary wedge transport time series
    wbw = dataset.createVariable('wbw',np.float64,(tdim.name))
    wbw.units = 'Sv'
    wbw.minimum_longitude = wbw_minlon
    wbw.maximum_longitude = int_minlon
    wbw.comment = 'Western boundary wedge transport time series (streamfunction at 1000m).'
    wbw[:] = wbw_trans.streamfunction[:,zind] / 1e6

    # Eastern-boundary wedge transport time series
    ebw = dataset.createVariable('ebw',np.float64,(tdim.name))
    ebw.units = 'Sv'
    ebw.minimum_longitude = int_maxlon
    ebw.maximum_longitude = ebw_maxlon
    ebw.comment = 'Eastern boundary wedge transport time series (streamfunction at 1000m).'
    ebw[:] = ebw_trans.streamfunction[:,zind] / 1e6

    # Upper mid ocean transport time series
    umo = dataset.createVariable('umo',np.float64,(tdim.name))
    umo.units = 'Sv'
    umo.minimum_longitude = wbw_minlon
    umo.maximum_longitude = ebw_maxlon
    umo.comment = 'Upper mid-ocean transport time series (streamfunction at 1000m). umo = wbw + geoint + ebw'
    umo[:] = wbw[:] + geoint[:] + ebw[:]

    # Meridional overturning transport time series - SAMBA approx
    moc_samba = dataset.createVariable('moc_samba',np.float64,(tdim.name))
    moc_samba.units = 'Sv'
    moc_samba.minimum_longitude = wbw_minlon
    moc_samba.maximum_longitude = ebw_maxlon
    moc_samba.comment = 'Time series of meridional overturning transport using SAMBA approximation (streamfunction at 1000m)'
    moc_samba[:] =  samba_trans.streamfunction[:,zind] / 1e6

    # Meridional overturning transport time series - model v
    moc_model = dataset.createVariable('moc_model',np.float64,(tdim.name))
    moc_model.units = 'Sv'
    moc_model.minimum_longitude = wbw_minlon
    moc_model.maximum_longitude = ebw_maxlon
    moc_model.comment = 'Time series of meridional overturning transport using model velocities (streamfunction at 1000m)'
    moc_model[:] =  model_trans.streamfunction[:,zind] / 1e6

    # Meridional overturning transport maxima time series - SAMBA approx
    mocmax_samba = dataset.createVariable('mocmax_samba',np.float64,(tdim.name))
    mocmax_samba.units = 'Sv'
    mocmax_samba.minimum_longitude = wbw_minlon
    mocmax_samba.maximum_longitude = ebw_maxlon
    mocmax_samba.comment = 'Time series of meridional overturning transport using SAMBA approximation (streamfunction maxima)'
    mocmax_samba[:] =  samba_trans.streamfunction.max(axis=1) / 1e6

    # Meridional overturning transport maxima time series - model v
    mocmax_model = dataset.createVariable('mocmax_model',np.float64,(tdim.name))
    mocmax_model.units = 'Sv'
    mocmax_model.minimum_longitude = wbw_minlon
    mocmax_model.maximum_longitude = ebw_maxlon
    mocmax_model.comment = 'Time series of meridional overturning transport using model velocities (streamfunction maxima)'
    mocmax_model[:] =  model_trans.streamfunction.max(axis=1) / 1e6

    # Overturning streamfunctions - SAMBA approx
    sf_samba = dataset.createVariable('sf_samba',np.float64,(tdim.name,zdim.name))
    sf_samba.units = 'Sv'
    sf_samba.minimum_longitude = wbw_minlon
    sf_samba.maximum_longitude = ebw_maxlon
    sf_samba.comment = 'Overturning streamfunctions using SAMBA approximation.'
    sf_samba[:] =  samba_trans.streamfunction/ 1e6

    # Meridional overturning transport time series - model v
    sf_model = dataset.createVariable('sf_model',np.float64,(tdim.name,zdim.name))
    sf_model.units = 'Sv'
    sf_model.minimum_longitude = wbw_minlon
    sf_model.maximum_longitude = ebw_maxlon
    sf_model.comment = 'Overturning streamfunctions using model velocities.'
    sf_model[:] =  model_trans.streamfunction / 1e6

    # Ekman stream function
    sf_ek = dataset.createVariable('sf_ek',np.float64,(tdim.name,zdim.name))
    sf_ek.units = 'Sv'
    sf_ek.minimum_longitude = int_minlon
    sf_ek.maximum_longitude = int_maxlon
    sf_ek.comment = 'Ekman overturning streamfunction.'
    sf_ek[:] =  ek_trans.streamfunction/ 1e6

    # wbw stream function
    sf_wbw = dataset.createVariable('sf_wbw',np.float64,(tdim.name,zdim.name))
    sf_wbw.units = 'Sv'
    sf_wbw.minimum_longitude = wbw_minlon
    sf_wbw.maximum_longitude = int_minlon
    sf_wbw.comment = 'Western boundary wedge overturning streamfunction.'
    sf_wbw[:] =  wbw_trans.streamfunction/ 1e6

    # ebw stream function
    sf_ebw = dataset.createVariable('sf_ebw',np.float64,(tdim.name,zdim.name))
    sf_ebw.units = 'Sv'
    sf_ebw.minimum_longitude = int_maxlon
    sf_ebw.maximum_longitude = ebw_maxlon
    sf_ebw.comment = 'Eastern boundary wedge overturning streamfunction.'
    sf_ebw[:] =  ebw_trans.streamfunction/ 1e6

    # Geostrophic interior stream function
    sf_geoint = dataset.createVariable('sf_geoint',np.float64,(tdim.name,zdim.name))
    sf_geoint.units = 'Sv'
    sf_geoint.minimum_longitude = int_minlon
    sf_geoint.maximum_longitude = int_maxlon
    sf_geoint.comment = 'Geostrophic interior overturning streamfunction.'
    sf_geoint[:] =  int_trans.streamfunction/ 1e6

    # model interior stream function
    sf_int_mod = dataset.createVariable('sf_int_mod',np.float64,(tdim.name,zdim.name))
    sf_int_mod.units = 'Sv'
    sf_int_mod.minimum_longitude = int_minlon
    sf_int_mod.maximum_longitude = int_maxlon
    sf_int_mod.comment = 'Model interior overturning streamfunction.'
    sf_int_mod[:] =  int_trans_mod.streamfunction/ 1e6

    # mid ocean stream function
    sf_mo = dataset.createVariable('sf_mo',np.float64,(tdim.name,zdim.name))
    sf_mo.units = 'Sv'
    sf_mo.minimum_longitude = int_minlon
    sf_mo.maximum_longitude = int_maxlon
    sf_mo.comment = 'Mid ocean overturning streamfunction (sf_mo = sf_wbw + sf_int + sf_ebw).'
    sf_mo[:] =  sf_geoint[:] + sf_wbw[:] + sf_ebw[:]

    # Total heat transport - SAMBA approx
    q_sum_samba = dataset.createVariable('q_sum_samba',np.float64,(tdim.name))
    q_sum_samba.units = 'PW'
    q_sum_samba.minimum_longitude = wbw_minlon
    q_sum_samba.maximum_longitude = ebw_maxlon
    q_sum_samba.comment = 'Total heat transport across section calculated using SAMBA approximations (q_sum_samba = q_ek + q_mo = q_ot_samba + q_gyre_samba + q_net_samba)'
    q_sum_samba[:] = samba_trans.oht_total / 1e15

    # Gyre heat transport - SAMBA approx
    q_gyre_samba = dataset.createVariable('q_gyre_samba',np.float64,(tdim.name))
    q_gyre_samba.units = 'PW'
    q_gyre_samba.minimum_longitude = wbw_minlon
    q_gyre_samba.maximum_longitude = ebw_maxlon
    q_gyre_samba.comment = 'Heat transport by the horizontal circulation calculated using SAMBA approximations '
    q_gyre_samba[:] = samba_trans.oht_by_horizontal / 1e15

    # Overturning heat transport - SAMBA approx
    q_ot_samba = dataset.createVariable('q_ot_samba',np.float64,(tdim.name))
    q_ot_samba.units = 'PW'
    q_ot_samba.minimum_longitude = wbw_minlon
    q_ot_samba.maximum_longitude = ebw_maxlon
    q_ot_samba.comment = 'Heat transport by the overturning circulation calculated using SAMBA approximations'
    q_ot_samba[:] = samba_trans.oht_by_overturning / 1e15

    # Heat transport by net throughflow - SAMBA approx
    q_net_samba = dataset.createVariable('q_net_samba',np.float64,(tdim.name))
    q_net_samba.units = 'PW'
    q_net_samba.minimum_longitude = wbw_minlon
    q_net_samba.maximum_longitude = ebw_maxlon
    q_net_samba.comment = 'Heat transport referenced to 0C by the net flow through the section using SAMBA approximations'
    q_net_samba[:] = samba_trans.oht_by_net / 1e15

    # Total heat transport - model v
    q_sum_model = dataset.createVariable('q_sum_model',np.float64,(tdim.name))
    q_sum_model.units = 'PW'
    q_sum_model.minimum_longitude = wbw_minlon
    q_sum_model.maximum_longitude = ebw_maxlon
    q_sum_model.comment = 'Total heat transport across section calculated using model velocities (q_sum_model = q_gyre_model + q_ot_model + q_net_model)'
    q_sum_model[:] = model_trans.oht_total / 1e15

    # Gyre heat transport -model v
    q_gyre_model = dataset.createVariable('q_gyre_model',np.float64,(tdim.name))
    q_gyre_model.units = 'PW'
    q_gyre_model.minimum_longitude = int_minlon
    q_gyre_model.maximum_longitude = int_maxlon
    q_gyre_model.comment = 'Heat transport by the horizontal circulation calculated using model velocities'
    q_gyre_model[:] = model_trans.oht_by_horizontal / 1e15

    # Overturning heat transport - model v
    q_ot_model = dataset.createVariable('q_ot_model',np.float64,(tdim.name))
    q_ot_model.units = 'PW'
    q_ot_model.minimum_longitude = wbw_minlon
    q_ot_model.maximum_longitude = ebw_maxlon
    q_ot_model.comment = 'Heat transport by the overturning circulation calculated using model velocities'
    q_ot_model[:] = model_trans.oht_by_overturning / 1e15

    # Heat transport by net throughflow - model v
    q_net_model = dataset.createVariable('q_net_model',np.float64,(tdim.name))
    q_net_model.units = 'PW'
    q_net_model.minimum_longitude = wbw_minlon
    q_net_model.maximum_longitude = ebw_maxlon
    q_net_model.comment = 'Heat transport referenced to 0C by the net flow through the section using model velocities'
    q_net_model[:] = model_trans.oht_by_net / 1e15

    # Heat transport by ekman
    q_ek = dataset.createVariable('q_ek',np.float64,(tdim.name))
    q_ek.units = 'PW'
    q_ek.minimum_longitude = int_minlon
    q_ek.maximum_longitude = int_maxlon
    q_ek.comment = 'Heat transport referenced to 0C by Ekman transport'
    q_ek[:] = ek_trans.oht_total / 1e15

    # Heat transport by wbw
    q_wbw = dataset.createVariable('q_wbw',np.float64,(tdim.name))
    q_wbw.units = 'PW'
    q_wbw.minimum_longitude = wbw_minlon
    q_wbw.maximum_longitude = int_minlon
    q_wbw.comment = 'Heat transport referenced to 0C by western boundary wedge transport'
    q_wbw[:] = wbw_trans.oht_total / 1e15

    # Heat transport by ebw
    q_ebw = dataset.createVariable('q_ebw',np.float64,(tdim.name))
    q_ebw.units = 'PW'
    q_ebw.minimum_longitude = int_maxlon
    q_ebw.maximum_longitude = ebw_maxlon
    q_ebw.comment = 'Heat transport referenced to 0C by eastern boundary wedge transport'
    q_ebw[:] = ebw_trans.oht_total / 1e15

    # Heat transport by zonal mean geostrophic interior
    q_geoint = dataset.createVariable('q_geoint',np.float64,(tdim.name))
    q_geoint.units = 'PW'
    q_geoint.minimum_longitude = int_minlon
    q_geoint.maximum_longitude = int_maxlon
    q_geoint.comment = 'Heat transport referenced to 0C by zonal mean of geostrophic interior transport'
    q_geoint[:] = (int_trans.oht_total - int_trans.oht_by_horizontal )/ 1e15

    # Heat transport by standing "eddy" component of geostrophic interior
    q_eddy = dataset.createVariable('q_eddy',np.float64,(tdim.name))
    q_eddy.units = 'PW'
    q_eddy.minimum_longitude = int_minlon
    q_eddy.maximum_longitude = int_maxlon
    q_eddy.comment = 'Heat transport referenced to 0C by standing eddy component of geostrophic interior transport'
    q_eddy[:] = (int_trans.oht_by_horizontal )/ 1e15

    # Heat transport by mid ocean
    q_mo = dataset.createVariable('q_mo',np.float64,(tdim.name))
    q_mo.units = 'PW'
    q_mo.minimum_longitude = wbw_minlon
    q_mo.maximum_longitude = ebw_maxlon
    q_mo.comment = 'Heat transport referenced to 0C by mid-ocean transport (q_mo = q_geoint + q_wbw + q_ebw + q_eddy)'
    q_mo[:] = q_geoint[:] + q_wbw[:] + q_ebw[:] + q_eddy[:]

    # Total freshwater transport - SAMBA approx
    fw_sum_samba = dataset.createVariable('fw_sum_samba',np.float64,(tdim.name))
    fw_sum_samba.units = 'Sv'
    fw_sum_samba.minimum_longitude = wbw_minlon
    fw_sum_samba.maximum_longitude = ebw_maxlon
    fw_sum_samba.reference_salinity = samba_trans.sref
    fw_sum_samba.comment = 'Total equivalent freshwater transport across section calculated using SAMBA approximations (fw_sum_samba = fw_fc + fw_ek + fw_mo = fw_ot_samba + fw_gyre_samba + fw_net_samba)'
    fw_sum_samba[:] = samba_trans.oft_total /1.0e6

    # Gyre freshwater transport - SAMBA approx
    fw_gyre_samba = dataset.createVariable('fw_gyre_samba',np.float64,(tdim.name))
    fw_gyre_samba.units = 'Sv'
    fw_gyre_samba.minimum_longitude = wbw_minlon
    fw_gyre_samba.maximum_longitude = ebw_maxlon
    fw_gyre_samba.comment = 'freshwater transport by the horizontal circulation calculated using SAMBA approximations '
    fw_gyre_samba[:] = samba_trans.oft_by_horizontal/1.0e6

    # Overturning freshwater transport - SAMBA approx
    fw_ot_samba = dataset.createVariable('fw_ot_samba',np.float64,(tdim.name))
    fw_ot_samba.units = 'Sv'
    fw_ot_samba.minimum_longitude = wbw_minlon
    fw_ot_samba.maximum_longitude = ebw_maxlon
    fw_ot_samba.comment = 'freshwater transport by the overturning circulation calculated using SAMBA approximations'
    fw_ot_samba[:] = samba_trans.oft_by_overturning /1.0e6

    # freshwater transport by net throughflow - SAMBA approx
    fw_net_samba = dataset.createVariable('fw_net_samba',np.float64,(tdim.name))
    fw_net_samba.units = 'Sv'
    fw_net_samba.minimum_longitude = wbw_minlon
    fw_net_samba.maximum_longitude = ebw_maxlon
    fw_net_samba.reference_salinity = samba_trans.sref
    fw_net_samba.comment = 'equivalent freshwater transport by the net flow through the section using SAMBA approximations'
    fw_net_samba[:] = samba_trans.oft_by_net /1.0e6

    # Total freshwater transport - model v
    fw_sum_model = dataset.createVariable('fw_sum_model',np.float64,(tdim.name))
    fw_sum_model.units = 'Sv'
    fw_sum_model.minimum_longitude = wbw_minlon
    fw_sum_model.maximum_longitude = ebw_maxlon
    fw_sum_model.reference_salinity = model_trans.sref
    fw_sum_model.comment = 'Total freshwater transport across section calculated using model velocities (fw_sum_model = fw_gyre_model + fw_ot_model + fw_net_model)'
    fw_sum_model[:] = model_trans.oft_total /1.0e6

    # Gyre freshwater transport -model v
    fw_gyre_model = dataset.createVariable('fw_gyre_model',np.float64,(tdim.name))
    fw_gyre_model.units = 'Sv'
    fw_gyre_model.minimum_longitude = wbw_minlon
    fw_gyre_model.maximum_longitude = ebw_maxlon
    fw_gyre_model.comment = 'freshwater transport by the horizontal circulation calculated using model velocities'
    fw_gyre_model[:] = model_trans.oft_by_horizontal /1.0e6

    # Overturning freshwater transport - model v
    fw_ot_model = dataset.createVariable('fw_ot_model',np.float64,(tdim.name))
    fw_ot_model.units = 'Sv'
    fw_ot_model.minimum_longitude = wbw_minlon
    fw_ot_model.maximum_longitude = ebw_maxlon
    fw_ot_model.comment = 'freshwater transport by the overturning circulation calculated using model velocities'
    fw_ot_model[:] = model_trans.oft_by_overturning /1.0e6

    # freshwater transport by net throughflow - model v
    fw_net_model = dataset.createVariable('fw_net_model',np.float64,(tdim.name))
    fw_net_model.units = 'Sv'
    fw_net_model.minimum_longitude = wbw_minlon
    fw_net_model.maximum_longitude = ebw_maxlon
    fw_net_model.reference_salinity = model_trans.sref
    fw_net_model.comment = 'equivalent freshwater transport by the net flow through the section using model velocities'
    fw_net_model[:] = model_trans.oft_by_net /1.0e6

    # freshwater transport by ekman
    fw_ek = dataset.createVariable('fw_ek',np.float64,(tdim.name))
    fw_ek.units = 'Sv'
    fw_ek.minimum_longitude = int_minlon
    fw_ek.maximum_longitude = int_maxlon
    fw_ek.reference_salinity = ek_trans.sref
    fw_ek.comment = 'equivalent freshwater transport by Ekman transport'
    fw_ek[:] = ek_trans.oft_total/1.0e6

    # freshwater transport by wbw
    fw_wbw = dataset.createVariable('fw_wbw',np.float64,(tdim.name))
    fw_wbw.units = 'Sv'
    fw_wbw.minimum_longitude = wbw_minlon
    fw_wbw.maximum_longitude = int_minlon
    fw_wbw.reference_salinity = wbw_trans.sref
    fw_wbw.comment = 'equivalent freshwater transport by western boundary wedge transport'
    fw_wbw[:] = wbw_trans.oft_total /1.0e6

    # freshwater transport by ebw
    fw_ebw = dataset.createVariable('fw_ebw',np.float64,(tdim.name))
    fw_ebw.units = 'Sv'
    fw_ebw.minimum_longitude = int_maxlon
    fw_ebw.maximum_longitude = ebw_maxlon
    fw_ebw.reference_salinity = ebw_trans.sref
    fw_ebw.comment = 'equivalent freshwater transport by eastern boundary wedge transport'
    fw_ebw[:] = ebw_trans.oft_total /1.0e6

    # freshwater transport by zonal mean geostrophic interior
    fw_geoint = dataset.createVariable('fw_geoint',np.float64,(tdim.name))
    fw_geoint.units = 'Sv'
    fw_geoint.minimum_longitude = int_minlon
    fw_geoint.maximum_longitude = int_maxlon
    fw_geoint.reference_salinity = int_trans.sref
    fw_geoint.comment = 'equivalent freshwater transport by zonal mean of geostrophic interior transport'
    fw_geoint[:] = (int_trans.oft_total - int_trans.oft_by_horizontal )/1.0e6

    # freshwater transport by standing "eddy" component of geostrophic interior
    fw_eddy = dataset.createVariable('fw_eddy',np.float64,(tdim.name))
    fw_eddy.units = 'Sv'
    fw_eddy.minimum_longitude = int_minlon
    fw_eddy.maximum_longitude = int_maxlon
    fw_eddy.reference_salinity = int_trans.sref
    fw_eddy.comment = 'equivalent freshwater transport by standing eddy component of geostrophic interior transport'
    fw_eddy[:] = (int_trans.oft_by_horizontal )/1.0e6

    # freshwater transport by mid ocean
    fw_mo = dataset.createVariable('fw_mo',np.float64,(tdim.name))
    fw_mo.units = 'Sv'
    fw_mo.minimum_longitude = int_minlon
    fw_mo.maximum_longitude = int_maxlon
    fw_mo.reference_salinity = int_trans.sref
    fw_mo.comment = 'equivalent freshwater transport by mid-ocean transport (fw_mo = fw_geoint + fw_wbw + fw_ebw + fw_eddy)'
    fw_mo[:] = fw_geoint[:] + fw_wbw[:] + fw_ebw[:] + fw_eddy[:]

    return dataset


