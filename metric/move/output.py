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
        suffix='_natl_meridional_transports_at_16N_td_endpoint.nc'
      else:
        suffix='_natl_meridional_transports_at_16N_td.nc'
    else:
      if config.getboolean('options', 'endpoint'):
        suffix='_natl_meridional_transports_at_16N_endpoint.nc'
      else:
        suffix='_natl_meridional_transports_at_16N.nc'

    savef = utils.get_savename(
        outdir, name, dates, date_format,
        suffix=suffix)
    dataset = Dataset(savef, 'w', format='NETCDF4_CLASSIC')

    return dataset

def create_netcdf(config, move_trans, model_trans, bdry_trans, int_trans, int_model_trans):
    """
    Return diagnostics for comparison with MOVE observations
    in netcdf object.

    """

    # Configuration options
    upper_level = config.getfloat('options','upper_level')
    lower_level = config.getfloat('options','lower_level')
    bdry_minlon = config.getfloat('options','bdry_minlon')
    bdry_maxlon = config.getfloat('options','bdry_maxlon')
    int_maxlon = config.getfloat('options','int_maxlon')
    georef = config.getfloat('options','georef_level')

    try:
        vref_level =  config.getfloat('options','vref_level')
    except:
        vref_level = 'None'

    # Create netcdf file and add dimensions
    dataset = open_ncfile(config, move_trans.dates)
    zdim = dataset.createDimension('z', move_trans.z.size)
    zbdim = dataset.createDimension('zbounds', move_trans.zbounds.size)
    xdim = dataset.createDimension('x', move_trans.x.size)
    xbdim = dataset.createDimension('xbounds', move_trans.xbounds.size)
    tdim = dataset.createDimension('time', None)

    # Create time variable
    time = dataset.createVariable('time',np.float64,(tdim.name,))
    time.units = 'hours since 0001-01-01 00:00:00.0'
    time.calendar = 'gregorian'
    time[:] = date2num(move_trans.dates, time.units, calendar=time.calendar)

    # Create z variable
    z = dataset.createVariable('z',np.float64,(zdim.name,))
    z.units = 'm'
    z[:] = move_trans.z

    # Create zbounds variable
    zbounds = dataset.createVariable('zbounds',np.float64,(zbdim.name,))
    zbounds.units = 'm'
    zbounds[:] = move_trans.zbounds

    # Create dz variable
    dz = dataset.createVariable('dz',np.float64,(zdim.name,))
    dz.units = 'm'
    dz[:] = move_trans.dz

    # Create x variable
    x = dataset.createVariable('x',np.float64,(xdim.name,))
    x.units = 'degrees_east'
    x[:] = move_trans.x

    # Create dx variable
    dx = dataset.createVariable('dx',np.float64,(xdim.name,))
    dx.units = 'm'
    dx[:] = move_trans.dx

    # Create xbounds variable
    xbounds = dataset.createVariable('xbounds',np.float64,(xbdim.name,))
    xbounds.units = 'degrees_east'
    xbounds[:] = move_trans.xbounds

    # Global attributes
    if config.getboolean('options', 'td_geo'):
      dataset.geostrophic_method = 'top-down'
    else:
      dataset.geostrophic_method = 'bottom-up'
      dataset.geostrophic_reference_level = georef
      dataset.reference_to_model_velocity = vref_level
    if config.getboolean('options', 'endpoint'):
      dataset.geostrophic_computation = 'endpoint'
    dataset.rhocp = move_trans.rhocp
    dataset.contact = 'fredc.ucar.edu'
    dataset.code_reference = 'https://github.com/NCAR/metric'
    dataset.method_references = 'Torsten Kanzow, Uwe Send, Walter Zenk, Alan D. Chave, Monika Rhein,\
 Monitoring the integrated deep meridional flow in the tropical North Atlantic: Long-term performance\
 of a geostrophic array, Deep Sea Research Part I: Oceanographic Research Papers, Volume 53, Issue 3,\
 2006, Pages 528-546, ISSN 0967-0637, https://doi.org/10.1016/j.dsr.2005.12.007.'

    # MOVE velocity cross section
    v = dataset.createVariable('v',np.float64,(tdim.name,zdim.name,xdim.name),fill_value=default_fillvals['f8'])
    v.units = 'm/s'
    v.comment = 'model velocity'
    v[:] = model_trans.v

    # MOVE geo velocity cross section
    vgeo = dataset.createVariable('vgeo',np.float64,(tdim.name,zdim.name,xdim.name),fill_value=default_fillvals['f8'])
    vgeo.units = 'm/s'
    vgeo.comment = 'geostrophic velocity'
    vgeo[:] = move_trans.v

    # MOVE dh cross section
    dh = dataset.createVariable('dh',np.float64,(tdim.name,zdim.name,xbdim.name),fill_value=default_fillvals['f8'])
    dh.units = 'm2/s2'
    dh.comment = 'dynamic heights'
    dh[:] = move_trans.dh

    # MOVE rho cross section
    rho = dataset.createVariable('rho',np.float64,(tdim.name,zdim.name,xbdim.name),fill_value=default_fillvals['f8'])
    rho.units = 'kg/m3'
    rho.comment = 'density'
    rho[:] = move_trans.rho

    # MOVE temp cross section
    temp = dataset.createVariable('temp',np.float64,(tdim.name,zdim.name,xdim.name),fill_value=default_fillvals['f8'])
    temp.units = 'degC'
    temp.comment = 'model potential temperature'
    temp[:] = model_trans.t

    # MOVE salinity cross section
    salt = dataset.createVariable('salt',np.float64,(tdim.name,zdim.name,xdim.name),fill_value=default_fillvals['f8'])
    salt.units = 'PSU'
    salt.comment = 'model salinity'
    salt[:] = model_trans.s

    # MOVE SSH cross section
    if config.getboolean('options', 'td_geo'):
      ssh = dataset.createVariable('ssh',np.float64,(tdim.name,xdim.name),fill_value=default_fillvals['f8'])
      ssh.units = 'm'
      ssh.comment = 'model SSH'
      ssh[:] = model_trans.ssh

    # MOVE potential temperature profile
    t_move = dataset.createVariable('t_move',np.float64,(tdim.name,zdim.name))
    t_move.units = 'degC'
    t_move.minimum_longitude = bdry_minlon
    t_move.maximum_longitude = int_maxlon
    t_move.comment = 'MOVE zonal mean potential temperature profile'
    t_move[:] = move_trans.zonal_avg_t

    # MOVE potential temperature profile
    t_int = dataset.createVariable('t_int',np.float64,(tdim.name,zdim.name))
    t_int.units = 'degC'
    t_int.minimum_longitude = bdry_maxlon
    t_int.maximum_longitude = int_maxlon
    t_int.comment = 'MOVE zonal mean potential temperature profile'
    t_int[:] = int_trans.zonal_avg_t

    # boundary component flow-weighted potential temperature
    t_bdry_fwt = dataset.createVariable('t_bdry_fwt',np.float64,(tdim.name))
    t_bdry_fwt.units = 'degC'
    t_bdry_fwt.minimum_longitude = bdry_minlon
    t_bdry_fwt.maximum_longitude = bdry_maxlon
    t_bdry_fwt.comment = 'boundary component flow-weighted potential temperature'
    t_bdry_fwt[:] = bdry_trans.oht_total / (bdry_trans.rhocp * bdry_trans.net_transport)

    # MOVE salinity profile
    s_move = dataset.createVariable('s_move',np.float64,(tdim.name,zdim.name))
    s_move.units = 'PSU'
    s_move.minimum_longitude = bdry_minlon
    s_move.maximum_longitude = int_maxlon
    s_move.comment = 'MOVE zonal mean salinity profile'
    s_move[:] = move_trans.zonal_avg_s

    # MOVE salinity profile
    s_int = dataset.createVariable('s_int',np.float64,(tdim.name,zdim.name))
    s_int.units = 'PSU'
    s_int.minimum_longitude = bdry_maxlon
    s_int.maximum_longitude = int_maxlon
    s_int.comment = 'MOVE zonal mean salinity profile'
    s_int[:] = int_trans.zonal_avg_s

    # boundary component flow-weighted salinity
    s_bdry_fwt = dataset.createVariable('s_bdry_fwt',np.float64,(tdim.name))
    s_bdry_fwt.units = 'PSU'
    s_bdry_fwt.minimum_longitude = bdry_minlon
    s_bdry_fwt.maximum_longitude = bdry_maxlon
    s_bdry_fwt.comment = 'boundary component flow-weighted salinity'
    s_bdry_fwt[:] = bdry_trans.oft_total / ((-1.0/bdry_trans.sref) * bdry_trans.net_transport)

    # MOVE transport profile - MOVE approx
    v_move = dataset.createVariable('v_move',np.float64,(tdim.name,zdim.name))
    v_move.units = 'Sv/m'
    v_move.minimum_longitude = bdry_minlon
    v_move.maximum_longitude = int_maxlon
    v_move.comment = 'MOVE transport profile using MOVE approximation'
    v_move[:] = move_trans.zonal_sum_v / 1e6

    # MOVE transport profile - model v
    v_model = dataset.createVariable('v_model',np.float64,(tdim.name,zdim.name))
    v_model.units = 'Sv/m'
    v_model.minimum_longitude = bdry_minlon
    v_model.maximum_longitude = int_maxlon
    v_model.comment = 'MOVE transport profile using model velocities'
    v_model[:] = model_trans.zonal_sum_v / 1e6

    # boundary component transport profile
    v_bdry = dataset.createVariable('v_bdry',np.float64,(tdim.name,zdim.name))
    v_bdry.units = 'Sv/m'
    v_bdry.minimum_longitude = bdry_minlon
    v_bdry.maximum_longitude = bdry_maxlon
    v_bdry.comment = 'boundary component transport profile'
    v_bdry[:] = bdry_trans.zonal_sum_v / 1e6

    # interior component transport profile
    v_int = dataset.createVariable('v_int',np.float64,(tdim.name,zdim.name))
    v_int.units = 'Sv/m'
    v_int.minimum_longitude = bdry_maxlon
    v_int.maximum_longitude = int_maxlon
    v_int.comment = 'interior component transport profile'
    v_int[:] = int_trans.zonal_sum_v / 1e6

    # interior component transport profile
    v_int_model = dataset.createVariable('v_int_model',np.float64,(tdim.name,zdim.name))
    v_int_model.units = 'Sv/m'
    v_int_model.minimum_longitude = bdry_maxlon
    v_int_model.maximum_longitude = int_maxlon
    v_int_model.comment = 'interior component transport profile using model velocities'
    v_int_model[:] = int_model_trans.zonal_sum_v / 1e6

    # Meridional overturning transport time series - MOVE approx
    moc_move = dataset.createVariable('moc_move',np.float64,(tdim.name))
    moc_move.units = 'Sv'
    moc_move.minimum_longitude = bdry_minlon
    moc_move.maximum_longitude = int_maxlon
    moc_move.comment = 'Time series of meridional overturning transport using MOVE approximation'
    moc_move[:] =  move_trans.net_transport / 1e6

    # Meridional overturning transport time series - model v
    moc_model = dataset.createVariable('moc_model',np.float64,(tdim.name))
    moc_model.units = 'Sv'
    moc_model.minimum_longitude = bdry_minlon
    moc_model.maximum_longitude = int_maxlon
    moc_model.comment = 'Time series of meridional overturning transport using model velocities'
    moc_model[:] =  model_trans.net_transport / 1e6

    # boundary component transport time series
    trans_bdry = dataset.createVariable('trans_bdry',np.float64,(tdim.name))
    trans_bdry.units = 'Sv'
    trans_bdry.minimum_longitude = bdry_minlon
    trans_bdry.maximum_longitude = bdry_maxlon
    trans_bdry.comment = 'boundary component transport time series'
    trans_bdry[:] = bdry_trans.net_transport / 1e6

    # internal component transport time series
    trans_int = dataset.createVariable('trans_int',np.float64,(tdim.name))
    trans_int.units = 'Sv'
    trans_int.minimum_longitude = bdry_maxlon
    trans_int.maximum_longitude = int_maxlon
    trans_int.comment = 'internal component transport time series'
    trans_int[:] = int_trans.net_transport / 1e6

    # internal component transport time series using model velocities
    trans_int_model = dataset.createVariable('trans_int_model',np.float64,(tdim.name))
    trans_int_model.units = 'Sv'
    trans_int_model.minimum_longitude = bdry_maxlon
    trans_int_model.maximum_longitude = int_maxlon
    trans_int_model.comment = 'internal component transport time series using model velocities'
    trans_int_model[:] = int_model_trans.net_transport / 1e6

    # Total heat transport - MOVE approx
    q_move = dataset.createVariable('q_move',np.float64,(tdim.name))
    q_move.units = 'PW'
    q_move.minimum_longitude = bdry_minlon
    q_move.maximum_longitude = int_maxlon
    q_move.comment = 'Total heat transport across section calculated using MOVE approximations'
    q_move[:] = move_trans.oht_total / 1e15

    # Total heat transport - model
    q_model = dataset.createVariable('q_model',np.float64,(tdim.name))
    q_model.units = 'PW'
    q_model.minimum_longitude = bdry_minlon
    q_model.maximum_longitude = int_maxlon
    q_model.comment = 'Total heat transport across section calculated using model velocities'
    q_model[:] = model_trans.oht_total / 1e15

    # Heat transport by boundary component
    q_bdry = dataset.createVariable('q_bdry',np.float64,(tdim.name))
    q_bdry.units = 'PW'
    q_bdry.minimum_longitude = bdry_minlon
    q_bdry.maximum_longitude = bdry_maxlon
    q_bdry.comment = 'Heat transport by the boundary component'
    q_bdry[:] = bdry_trans.oht_total / 1e15

    # Heat transport by internal component
    q_int = dataset.createVariable('q_int',np.float64,(tdim.name))
    q_int.units = 'PW'
    q_int.minimum_longitude = bdry_maxlon
    q_int.maximum_longitude = int_maxlon
    q_int.comment = 'Heat transport by the internal component'
    q_int[:] = int_trans.oht_total / 1e15

    # internal component heat transport -model v
    q_int_model = dataset.createVariable('q_int_model',np.float64,(tdim.name))
    q_int_model.units = 'PW'
    q_int_model.minimum_longitude = bdry_maxlon
    q_int_model.maximum_longitude = int_maxlon
    q_int_model.comment = 'Heat transport by the internal component using model velocities'
    q_int_model[:] = (model_trans.oht_total - bdry_trans.oht_total) / 1e15

    # Total freshwater transport - MOVE approx
    fw_move = dataset.createVariable('fw_move',np.float64,(tdim.name))
    fw_move.units = 'Sv'
    fw_move.minimum_longitude = bdry_minlon
    fw_move.maximum_longitude = int_maxlon
    fw_move.reference_salinity = move_trans.sref
    fw_move.comment = 'Total equivalent freshwater transport across section calculated using MOVE approximations'
    fw_move[:] = move_trans.oft_total /1.0e6

    # Total freshwater transport - model v
    fw_model = dataset.createVariable('fw_model',np.float64,(tdim.name))
    fw_model.units = 'Sv'
    fw_model.minimum_longitude = bdry_minlon
    fw_model.maximum_longitude = int_maxlon
    fw_model.reference_salinity = model_trans.sref
    fw_model.comment = 'Total freshwater transport across section calculated using model velocities'
    fw_model[:] = model_trans.oft_total /1.0e6

    # freshwater transport by boundary component
    fw_bdry = dataset.createVariable('fw_bdry',np.float64,(tdim.name))
    fw_bdry.units = 'Sv'
    fw_bdry.minimum_longitude = bdry_minlon
    fw_bdry.maximum_longitude = bdry_maxlon
    fw_bdry.reference_salinity = bdry_trans.sref
    fw_bdry.comment = 'freshwater transport by the boundary component'
    fw_bdry[:] = bdry_trans.oft_total /1.0e6

    # internal component freshwater transport - MOVE approx
    fw_int = dataset.createVariable('fw_int',np.float64,(tdim.name))
    fw_int.units = 'Sv'
    fw_int.minimum_longitude = bdry_maxlon
    fw_int.maximum_longitude = int_maxlon
    fw_int.comment = 'freshwater transport by the internal component'
    fw_int[:] = int_trans.oft_by_horizontal/1.0e6

    # internal component freshwater transport -model v
    fw_int_model = dataset.createVariable('fw_int_model',np.float64,(tdim.name))
    fw_int_model.units = 'Sv'
    fw_int_model.minimum_longitude = bdry_maxlon
    fw_int_model.maximum_longitude = int_maxlon
    fw_int_model.comment = 'freshwater transport by the internal component using model velocities'
    fw_int_model[:] = (model_trans.oft_by_horizontal - bdry_trans.oft_total) /1.0e6

    return dataset


