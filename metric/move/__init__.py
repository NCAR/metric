"""
METRIC is a package to calculate diagnostics of the Atlantic Meridional 
Overturning Circulation (AMOC) using output from an ocean general circulation
model for comparison with observed data from the MOVE array at 16N.

METRIC is adapted from RapidMoc by Chris Roberts at ECMWF (chris.roberts@ecmwf.int)

Author:
Fred Castruccio, NCAR (fredc@ucar.edu)
"""

from .transports import *
from .output import *
from .plotdiag import *
from .observations import *

