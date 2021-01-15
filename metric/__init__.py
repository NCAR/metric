"""
METRIC is a package to calculate diagnostics of the Atlantic Meridional 
Overturning Circulation (AMOC) using output from an ocean general circulation
model for comparison with observed data from the MOVE array at 16N.

METRIC is adapted from RapidMoc by Chris Roberts at ECMWF (chris.roberts@ecmwf.int)

Author:
Fred Castruccio, NCAR (fredc@ucar.edu)
"""


from pkg_resources import DistributionNotFound, get_distribution

from .metric import *
#from .utils import *
#from .sections import *
#from .geostrophy import *
import metric.utils
import metric.sections
import metric.geostrophy

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass
