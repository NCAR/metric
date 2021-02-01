"""
Module containing code to plot transport diagnostics for SAMBA

"""


import matplotlib
#matplotlib.use('AGG')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm

import numpy as np
from scipy import stats

from metric.samba import observations
from metric import utils


# COLORS 
c1='#a6cee3'
c2='#1f78b4'
c3='#b2df8a'
c4='#33a02c'
MAX_NAME_LEN=35


def plot_volume_components(config, trans, basename='', name='simulated', obs=None, lw=4):
    """ Plot volume transport component time series """
    
    # Add model data to sub-axis
    fig = plt.figure(figsize=(8,12))
    fig.add_subplot(2,1,1)

    dts = utils.get_ncdates(config, trans)
    wbw = trans.variables['wbw'][:]
    ebw = trans.variables['ebw'][:]
    geoint = trans.variables['geoint'][:]
    ekman = trans.variables['ekman'][:]
    moc = trans.variables['moc_samba'][:]

    wbw_label = 'WBW transport (%6.1f Sv)' % (wbw.mean())
    ebw_label = 'EBW transport (%6.1f Sv)' % (ebw.mean())
    geoint_label = 'Geo. int. transport (%6.1f Sv)' % (geoint.mean())
    ekman_label = 'Ekman transport (%6.1f Sv)' % (ekman.mean())
    moc_label = 'Total MOC SAMBA (%6.1f Sv)' % (moc.mean())

    plt.plot(dts, wbw, linewidth=lw, color=c1, label=wbw_label)
    plt.plot(dts, ebw, linewidth=lw, color=c2, label=ebw_label)
    plt.plot(dts, geoint, linewidth=lw, color=c3, label=geoint_label)
    plt.plot(dts, ekman, linewidth=lw, color=c4, label=ekman_label)
    plt.plot(dts, moc, linewidth=lw, color='k', label=moc_label)

    plt.xlabel('Date')
    plt.ylim([-50,50])
    plt.ylabel('Sverdrups')
    plt.title('Overturning components at 34.5S in %s' % name)
    plt.legend(loc=8, fontsize=8, ncol=2)


    # Add optional observational data to sub-axis
    if (obs is not None):
        fig.add_subplot(2,1,2)
        
        moc_obs_label = 'MOC (%6.1f Sv)' % (obs.moc[~np.isnan(obs.moc)].mean())
        
        plt.plot(obs.dates, obs.moc, linewidth=lw, color='k', label=moc_obs_label)

        plt.xlabel('Date')
        plt.ylim([-50,50])
        plt.ylabel('Sverdrups')
        plt.title('Overturning components at 34.5S in SAMBA observations')
        plt.legend(loc=8, fontsize=8, ncol=2)

    # Save plot
    plt.tight_layout()
    savef = basename + 'volume_transport_components_at_34.5s.png'
    print('SAVING: {}'.format(savef))
    fig.savefig(savef, dpi=300)
    plt.close()
    
        

def plot_diagnostics(config, trans):
    """ Plot volume and heat transport diagnostics against SAMBA observations """

    # Initialize data    
    name='simulated'
    outdir='./'
    date_format='%Y%m%d',
    lw=4

    # Initialize observations
    obs_trans_f = None

    # Get observation file paths
    time_avg = utils.get_config_opt(config, 'observations', 'time_avg')
    obs_trans_f = utils.get_config_opt(config, 'observations', 'volume_transports')

    # Load observations, if specified
    if obs_trans_f is not None:
        obs_trans = observations.TransportObs(obs_trans_f, time_avg=time_avg)

    # Call plot routines
    outdir = config.get('output', 'outdir')
    date_format = config.get('output', 'date_format')
    name = config.get('output', 'name')

    # Create basename for output files
    dts = utils.get_ncdates(config, trans)
    basename = utils.get_savename(outdir, name, dts, date_format,suffix='_')
    name = name[0:MAX_NAME_LEN]

    # Plot data
    plot_volume_components(config, trans, basename=basename, name=name, obs=obs_trans, lw=lw)
   
