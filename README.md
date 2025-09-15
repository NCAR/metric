# **Meridional ovErTurning ciRculation diagnostIC (METRIC)**

METRIC is a fork of the [RapidMoc](https://github.com/cdr30/RapidMoc) package, which extends the calculation of observation-style transports to other observing arrays.

The METRIC python package enables consistent calculations of Atlantic meridional overturning circulation (AMOC) 
estimates at various observational arrays from ocean general circulation models. To make the most appropriate comparisons, the package evaluates the model meridional overturning circulation using analogous observation-style methods.

The current version allows AMOC estimates at the RAPID (26.5N) array, the MOVE (16N) array, and the SAMBA (34.5S) array. METRIC also includes a few additional, alternative approaches to calculate these transports.

In this latest version of METRIC, the package has been updated to utilise xarray and dask to allow users to calculate AMOC estimates from large multi-file datasets.

Installing METRIC
-------------

**Downloading the code via SSH:**
```shell
   git clone git@github.com:NCAR/metric.git
```

**Installing the updated METRIC package:**

```shell
   cd metric
   pip install -e .
```
Using METRIC
------------
**Running the code:**
```shell
   ./run_metric.py [-c] [-t] [-s] [-v] [-ssh] [-taux]
```

   Calculate AMOC estimates at the RAPID, MOVE, and SAMBA arrays using the observation-style method.

   ### Positional Arguments:
   **-c** = Path to configuration file.\
   **-t** = Path to netcdf file containing temperature data.\
   **-s** = Path to netcdf file containing salinity data.\
   **-v** = Path to netcdf file containing meridional velocity data.\
   **-ssh** = Path to netcdf file containing sea surface height data.\
   **-taux** = Path to netcdf file containing zonal wind stress data.

   ### Optional Arguments:
   **--outdir** = Path to output directory (overwrites configuration file).\
   **--name** = Name of the output file (overwrites configuration file).\
   **--shift** = Shift output dates for plotting.

Citation
--------

Please cite the associated digital object identifier and Danabasoglu et al. (2021) if you use this module for your research.

Castruccio F. S., 2021: NCAR/metric: metric v0.1. doi/10.5281/zenodo.4708277

[![DOI](https://zenodo.org/badge/331704850.svg)](https://zenodo.org/badge/latestdoi/331704850)

Danabasoglu, G., Castruccio, F. S.,  Small, R. J., Tomas, R., Frajka-Williams, E., and Lankhorst, M., (2021). Revisiting AMOC Transport Estimates from Observations and Models. *Geophysical Research Letters*, in review. 

The original RapidMoc package should be cited using the associated digital object identifiers and [Roberts et al. (2013)](http://onlinelibrary.wiley.com/doi/10.1002/grl.50930/full)

Roberts, C.D., 2017: cdr30/RapidMoc: RapidMoc v1.0.1. doi:10.5281/zenodo.1036387

Roberts, C. D., et al. (2013), Atmosphere drives recent interannual variability of the Atlantic meridional overturning circulation at 26.5°N, Geophys. Res. Lett., 40, 5164–5170 doi:10.1002/grl.50930.
