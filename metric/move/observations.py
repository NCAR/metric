"""
Module containing code to work with Move observational data

"""

from netCDF4 import Dataset, num2date, date2num
import datetime
import numpy as np
import metric.utils

class MoveObs(object):
    """ Template class to interface with observed ocean transports """

    def __init__(self, f, time_avg=None, mindt=None, maxdt=None):
        """ Create instance holding ocean transport data """
        self.f = f
        self.time_avg = time_avg
        self.mindt = mindt
        self.maxdt = maxdt
        self._read_data()

    def _read_data(self):
        """ Abstract method to read data and apply time averaging """
        pass
    
    def _read_dates(self):
        """ Abstract method to initialized dates """
        pass
    
    def _ym_dates(self):
        """ Return yearly mean date time objects """
        ym_dates = []
        
        for yr in range(self.yy.min(), self.yy.max()+1):
            ind = (self.yy == yr)
            if ind.any():
                ym_dates.append(datetime.datetime(yr, 7, 1))
                    
        return np.array(ym_dates)
    
    def _mm_dates(self):
        """ Return monthly mean date time objects """
        mm_dates = []
        
        for yr in range(self.yy.min(), self.yy.max()+1):
            for mon in range(1,12+1):
                ind = (self.yy == yr) & (self.mm == mon)
                if ind.any():
                    mm_dates.append(datetime.datetime(yr, mon, 15))
                    
        return np.array(mm_dates)
    
    def _calc_ym(self, data, profile=False):
        """ Return yearly mean values """
        ym_data = []
        
        for yr in range(self.yy.min(), self.yy.max()+1):
            ind = (self.yy == yr)
            if ind.any():
                if profile:
                    ym_data.append(np.mean(data[ind,:],axis=0))
                else:
                    ym_data.append(np.mean(data[ind]))
                    
        return np.array(ym_data)
    
    def _calc_mm(self, data, profile=False, profile_bdry=False):
        """ Return monthly mean values """
        mm_data = []
        
        for yr in range(self.yy.min(), self.yy.max()+1):
            for mon in range(1,12+1):
                ind = (self.yy == yr) & (self.mm == mon)
                if ind.any():
                    if profile:
                        mm_data.append(np.mean(data[ind,:],axis=0))
                    elif profile_bdry:
                        mm_data.append(np.mean(data[:,ind,:],axis=1))
                    else:
                        mm_data.append(np.mean(data[ind]))
        
        return np.array(mm_data)
         
    def _readnc(self, ncvar):
        """ Read variable from netcdf file """
        nc = Dataset(self.f)
        data = nc.variables[ncvar][:]
        nc.close()
        
        return data
    

class TransportObs(MoveObs):
    """ 
    Sub-class to hold volume transport observations
    from the MOVE array at 16N.
    
    Data source:
    http://mooring.ucsd.edu/dev/move/
    
    Data reference:
    http://dx.doi.org/10.1016/j.dsr.2005.12.007 
    http://dx.doi.org/10.1029/2011GL049801
 
    """
     
    def _read_data(self):
        """ Read data and apply time averaging """
        self._read_dates()
        
        if self.time_avg is None:
            self.dates = self.original_dates
            self.trans_total = self._readnc('TRANSPORT_TOTAL')
            self.trans_int = self._readnc('transport_component_internal')
            self.trans_int_offset = self._readnc('transport_component_internal_offset')
            self.trans_bdry = self._readnc('transport_component_boundary')           
        elif self.time_avg == 'monthly':
            self.dates = self._mm_dates()
            self.trans_total = self._calc_mm(self._readnc('TRANSPORT_TOTAL'))
            self.trans_int = self._calc_mm(self._readnc('transport_component_internal'))
            self.trans_int_offset = self._calc_mm(self._readnc('transport_component_internal_offset'))
            self.trans_bdry = self._calc_mm(self._readnc('transport_component_boundary'))           
        elif self.time_avg == 'yearly':
            self.dates = self._ym_dates()
            self.trans_total = self._calc_ym(self._readnc('TRANSPORT_TOTAL'))
            self.trans_int = self._calc_ym(self._readnc('transport_component_internal'))
            self.trans_int_offset = self._calc_ym(self._readnc('transport_component_internal_offset'))
            self.trans_bdry = self._calc_ym(self._readnc('transport_component_boundary'))           
        else:
            print(self.time_avg)
            raise ValueError('time_avg must be "monthly" or "yearly"')
        
        if (self.mindt is not None) and (self.maxdt is not None):
            tind = utils.get_dateind(self.dates, self.mindt, self.maxdt)
            self.trans_total = self.trans_total[tind]
            self.trans_int = self.trans_int[tind]
            self.trans_int_offset = self.trans_int_offset[tind]
            self.trans_bdry = self.trans_bdry[tind]
            self.dates = self.dates[tind]

    def _read_dates(self):
        """ Read date information from file """
        nc = Dataset(self.f)
        t = nc.variables['TIME']
        self.original_dates = num2date(t[:],units=t.units)
        self.hh = np.array([dt.hour for dt in self.original_dates], dtype=np.int)
        self.dd = np.array([dt.day for dt in self.original_dates], dtype=np.int)
        self.mm = np.array([dt.month for dt in self.original_dates], dtype=np.int)
        self.yy = np.array([dt.year for dt in self.original_dates], dtype=np.int)


