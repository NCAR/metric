"""
Module containing code to work with SAMBA observational data

"""

from netCDF4 import Dataset, num2date, date2num
import datetime
import numpy as np
import metric.utils

class SambaObs(object):
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

    def _calc_mm(self, data, profile=False):
        """ Return monthly mean values """
        mm_data = []

        for yr in range(self.yy.min(), self.yy.max()+1):
            for mon in range(1,12+1):
                ind = (self.yy == yr) & (self.mm == mon)
                if ind.any():
                    if profile:
                        mm_data.append(np.mean(data[ind,:],axis=0))
                    else:
                        mm_data.append(np.mean(data[ind]))

        return np.array(mm_data)

    def _readnc(self, ncvar):
        """ Read variable from netcdf file """
        nc = Dataset(self.f)
        data = nc.variables[ncvar][:]
        nc.close()

        return data



class TransportObs(SambaObs):
    """
    Sub-class to hold volume transport observations
    from the SAMBA array at 34.5S.

    Data source:
    https://www.aoml.noaa.gov/phod/SAMOC_international/samoc_data.php

    Data reference:
    Meinen, C. S., Speich, S., Piola, A. R., Ansorge, I., Campos, E., Kersale, M., et al. (2018).
    Meridional overturning circula- tion transport variability at 34.5°S during 200-2017:
    Baroclinic and barotropic flows and the dueling influence of the boundaries.
    Geophysical Research Letters, 45, 4180-4188. https://doi.org/ 10.1029/2018GL077408

    """

    def _read_data(self):
        """ Read data and apply time averaging """
        self._read_dates()

        if self.time_avg is None:
            self.dates = self.original_dates
            self.moc_clinic = self._readnc('moc_clinic')
            self.moc_tropic = self._readnc('moc_tropic')
            self.moc_ek = self._readnc('moc_ek')
            self.moc_west_dens = self._readnc('moc_west_dens')
            self.moc_east_dens = self._readnc('moc_east_dens')
            self.moc_west_pres = self._readnc('moc_west_pres')
            self.moc_east_pres = self._readnc('moc_east_pres')
            self.moc = self._readnc('moc')
        elif self.time_avg == 'monthly':
            self.dates = self._mm_dates()
            self.moc_clinic = self._calc_mm(self._readnc('moc_clinic'))
            self.moc_tropic = self._calc_mm(self._readnc('moc_tropic'))
            self.moc_ek = self._calc_mm(self._readnc('moc_ek'))
            self.moc_west_dens = self._calc_mm(self._readnc('moc_west_dens'))
            self.moc_east_dens = self._calc_mm(self._readnc('moc_east_dens'))
            self.moc_west_pres = self._calc_mm(self._readnc('moc_west_pres'))
            self.moc_east_pres = self._calc_mm(self._readnc('moc_east_pres'))
            self.moc = self._calc_mm(self._readnc('moc'))
        elif self.time_avg == 'yearly':
            self.dates = self._ym_dates()
            self.moc_clinic = self._calc_ym(self._readnc('moc_clinic'))
            self.moc_tropic = self._calc_ym(self._readnc('moc_tropic'))
            self.moc_ek = self._calc_ym(self._readnc('moc_ek'))
            self.moc_west_dens = self._calc_ym(self._readnc('moc_west_dens'))
            self.moc_east_dens = self._calc_ym(self._readnc('moc_east_dens'))
            self.moc_west_pres = self._calc_ym(self._readnc('moc_west_pres'))
            self.moc_east_pres = self._calc_ym(self._readnc('moc_east_pres'))
            self.moc = self._calc_ym(self._readnc('moc'))
        else:
            print(self.time_avg)
            raise ValueError('time_avg must be "monthly" or "yearly"')

        if (self.mindt is not None) and (self.maxdt is not None):
            tind = utils.get_dateind(self.dates, self.mindt, self.maxdt)
            self.moc_clinic = self.moc_clinic[tind]
            self.moc_tropic = self.moc_tropic[tind]
            self.moc_ek = self.moc_ek[tind]
            self.moc_west_dens = self.moc_west_dens[tind]
            self.moc_east_dens = self.moc_east_dens[tind]
            self.moc_west_pres = self.moc_west_pres[tind]
            self.moc_east_pres = self.moc_east_pres[tind]
            self.moc = self.moc[tind]

    def _read_dates(self):
        """ Read date information from file """
        nc = Dataset(self.f)
        t = nc.variables['time']
        self.original_dates = num2date(t[:],units=t.units)
        self.hh = np.array([dt.hour for dt in self.original_dates], dtype=np.int)
        self.dd = np.array([dt.day for dt in self.original_dates], dtype=np.int)
        self.mm = np.array([dt.month for dt in self.original_dates], dtype=np.int)
        self.yy = np.array([dt.year for dt in self.original_dates], dtype=np.int)


