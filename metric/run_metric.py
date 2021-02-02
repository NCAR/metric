#!/usr/bin/env python

"""
metric command line script
"""


import sys
from .metric import compute_amoc_transport



def main():

  args = dict(zip(sys.argv[1::2],sys.argv[2::2]))

  print('')
  print('')
  if 'move' in args['-c']:
    print('Run MOVE AMOC transport computation using:')
  elif 'rapid' in args['-c']:
    print('Run RAPID AMOC transport computation using:')
  elif 'samba' in args['-c']:
    print('Run SAMBA SAMOC transport computation using:')
  for key in args.keys(): 
    if key == '-c':
      print('Path to config file: {}'.format(args[key]))
    if key == '-t':
      print('Path to temperature file: {}'.format(args[key]))
    if key == '-s':
      print('Path to salinity file: {}'.format(args[key]))
    if key == '-v':
      print('Path to meridional velocity file: {}'.format(args[key]))
    if key == '-ssh':
      print('Path to SSH file: {}'.format(args[key]))
    if key == '-taux':
      print('Path to zonal wind stress file: {}'.format(args[key]))
    if key == '--name':
      print('Output files name (Overrides value in config file): {}'.format(args[key]))
    if key == '--outdir':
      print('Output data path (Overrides value in config file): {}'.format(args[key]))
    if key == '--shift':
      print('Shift output dates for plotting such that the output time series start with {}'.format(args[key]))

  compute_amoc_transport(sys.argv[1:])

  print('')
  print('Done')

