"""
Module containing utility functions

"""


import argparse as ap
import configparser
import copy
from os import path
import errno

import numpy as np
import datetime
import cftime


class ShapeError(Exception):
    pass


def get_args(line_args=None):
    """   Get arguments from command line.  """

    parser = ap.ArgumentParser(
        description='Gathers from command line to calculate AMOC diagnostics using ocean model data.')
    parser.add_argument('-c', type=str, action='store', dest='config_file',
                        required=True, help='Path for netcdf file(s) containing temperature data.')
    parser.add_argument('-t', type=str, action='store', dest='temp_file',
                        required=True, help='Path for netcdf file(s) containing temperature data.')
    parser.add_argument('-s', type=str, action='store', dest='salt_file',
                        required=True, help='Path for netcdf file(s) containing salinity data.')
    parser.add_argument('-v', type=str, action='store', dest='vel_file',
                        required=True, help='Path for netcdf file(s) containing meridional velocity data.')
    parser.add_argument('-taux', type=str, action='store', dest='taux_file',
                        required=False, help='Path for netcdf file(s) containing zonal wind stress data.')
    parser.add_argument('-ssh', type=str, action='store', dest='ssh_file',
                        required=False, help='Path for netcdf file(s) containing Sea Surface Height data.')
    parser.add_argument('--name', type=str, action='store', dest='name', 
                        default=None, help='Name used in output files. Overrides value in config file.')
    parser.add_argument('--outdir', type=str, action='store', dest='outdir',
                        default=None, help='Path used for output data. Overrides value in config file.')
    parser.add_argument('--shift', type=str, action='store', dest='shift_date',
                        default=None, help='Shift dates.')
    # Gather the provided arguements as an array.
    args = parser.parse_args(line_args)

    return args


def get_config(args):
    """ Return configuration options as <configparser> object. """
    config = configparser.ConfigParser()
    config.read(args.config_file)

    return config


def get_config_opt(config, section, option):
    """ Return option if exists, else None """
    if config.has_option(section, option):
        return config.get(section, option)
    else:
        return None


def get_ncdates(config, nc, tvar='time'):
    """ Return dates from netcdf time coordinate """
    t = nc.variables[tvar]
    dts = cftime.num2date(t[:], t.units, calendar=t.calendar)
    mpldts = np.array([datetime.datetime(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second) for dt in dts])
    if config.has_option('output', 'shift_date'):
      dtshift = config.get('output', 'shift_date')
      mpldtshift = datetime.datetime(np.int(dtshift[:4]), np.int(dtshift[4:6]), np.int(dtshift[6:8]))
      tdelta = mpldtshift - mpldts[0]
      mpldts = mpldts + tdelta
      
    return mpldts


def get_datestr(dates, date_format):
    """ Return string of form 'mindate-maxdate' for specified format """
    datestr = '%s-%s' % (dates[0].strftime(date_format).replace(' ', '0'), 
                         dates[-1].strftime(date_format).replace(' ', '0'))
    
    return datestr


def get_savename(outdir, name, dates, date_format, suffix=''):
    """ Return savename for output file """
    datestr = get_datestr(dates, date_format)
    savename = path.join(outdir,'%s_%s%s' % (name, datestr, suffix))
    
    return savename


def get_daterange(dates1, dates2):
    """ Return min and max date range for overlapping period """
    mindt = max([dates1.min(),dates2.min()])
    maxdt = min([dates1.max(),dates2.max()])
    
    return mindt, maxdt


def get_dateind(dates, mindt, maxdt):
    """ Return index to extract data over specified date range """
    return (dates >= mindt) & (dates <= maxdt)


def get_indrange(vals,minval,maxval):
    """ 
    Return max and min indices to use for simple slicing when
    updating multi-dim np.array that avoids creation of copies.
        
    """ 
    
    if vals.ndim != 1:
        raise ShapeError('get_inds: expected 1-d numpy array')
        
    inds = np.where((vals >= minval) & (vals < maxval))[0]
    
    if len(inds) != 0:
        minind = inds.min()
        maxind = inds.max() + 1
    else:
        raise ValueError('get_inds: no matching data')
        
    return minind, maxind
    
    

def find_nearest(array, val, min=False):
    """ Returns index for value in array that is closest to val. """
    if not isinstance(array, np.ndarray):
        raise TypeError('Array must be a np.ndarray object')

    abs_array = np.absolute(array - val)
    minval = abs_array.min()
    inds = np.where(abs_array == minval)[0]

    if len(inds) > 1:
        if min:
            inds = np.where(array == array[inds].min())[0]
        else:
            inds = np.where(array == array[inds].max())[0]
    
    return inds
