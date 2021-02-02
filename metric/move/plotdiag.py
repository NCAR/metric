"""
Module containing code to plot transport diagnostics for MOVE

"""


import matplotlib
#matplotlib.use('AGG')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm

import numpy as np
from scipy import stats

from metric.move import observations
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
    trans_bdry = trans.variables['trans_bdry'][:]
    trans_int = trans.variables['trans_int'][:]
    moc = trans.variables['moc_move'][:]

    bdry_label = 'Boundary transport (%6.1f Sv)' % (trans_bdry.mean())
    int_label = 'Internal transport (%6.1f Sv)' % (trans_int.mean())
    moc_label = 'MOC (%6.1f Sv)' % (moc.mean())

    plt.plot(dts, trans_bdry, linewidth=lw, color=c1, label=bdry_label)
    plt.plot(dts, trans_int, linewidth=lw, color=c2, label=int_label)
    plt.plot(dts, moc, linewidth=lw, color='k', label=moc_label)

    plt.xlabel('Date')
    plt.ylim([-50,50])
    plt.ylabel('Sverdrups')
    plt.title('Overturning components at 16N in %s' % name)
    plt.legend(loc=8, fontsize=8, ncol=2)


    # Add optional observational data to sub-axis
    if (obs is not None):
        fig.add_subplot(2,1,2)
        
        trans_bdry_obs_label = 'Boundary transport (%6.1f Sv)' % (obs.trans_bdry[~np.isnan(obs.trans_bdry)].mean())
        trans_int_obs_label = 'Internal transport (%6.1f Sv)' % (obs.trans_int[~np.isnan(obs.trans_int)].mean())
        trans_int_offset_obs_label = 'Internal transport offset (%6.1f Sv)' % (obs.trans_int_offset[~np.isnan(obs.trans_int_offset)].mean())
        trans_total_obs_label = 'MOC (%6.1f Sv)' % (obs.trans_total[~np.isnan(obs.trans_total)].mean())
        
        plt.plot(obs.dates, obs.trans_bdry, linewidth=lw, color=c1, label=trans_bdry_obs_label)
        plt.plot(obs.dates, obs.trans_int, linewidth=lw, color=c2, label=trans_int_obs_label)
        plt.plot(obs.dates, obs.trans_int_offset, linewidth=lw, color=c3, label=trans_int_offset_obs_label)
        plt.plot(obs.dates, obs.trans_total, linewidth=lw, color='k', label=trans_total_obs_label)

        plt.xlabel('Date')
        plt.ylim([-50,50])
        plt.ylabel('Sverdrups')
        plt.title('Overturning components at 16N in MOVE observations')
        plt.legend(loc=8, fontsize=8, ncol=2)

    # Save plot
    plt.tight_layout()
    savef = basename + 'volume_transport_components_at_16n.png'
    print('SAVING: {}'.format(savef))
    fig.savefig(savef, dpi=300)
    plt.close()
    
        

def plot_diagnostics(config, trans):
    """ Plot volume and heat transport diagnostics against MOVE observations """

    # Initialize data    
    name='simulated'
    outdir='./'
    date_format='%Y%m%d',
    lw=4

    # Initialize observations
    obs_trans_f, obs_vel_f = None, None

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
   
