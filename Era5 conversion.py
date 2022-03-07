"""
ERA5 reanalysis data is at hourly timescale; ; represented in 'm of water equivalent per hour'
(m of water eq/hr)
"""
E_hourly = (xr.open_dataset('/home/chandra/data/Kruger/era5/evap_ERA5_reanalysis_hourly_2007.nc').e[:,::-1,:])
E_hourly_to_monthly = ((E_hourly.resample(time = '1M').sum('time'))*(-1000))
E_hourly_to_year = ((E_hourly_to_monthly.sum('time')))
#E_hourly_to_year.plot()
E_hourly_to_monthly.to_netcdf('/home/chandra/data/Kruger/era5/evap_ERA5_reanalysis_hourly_2007-resampled_monthly.nc')

"""
ERA5 monthly averaged reanalysis data is at montly timescale; ; represnted in 'm of water equivalent per day'
(m of water eq/day)
"""
E_avgmonthly = (xr.open_dataset('/home/chandra/data/Kruger/era5/evap_ERA5_monthly_averaged_reanalysis_2007.nc').e[:,::-1,:])
days_in_month = [31,28,31,30,31,30,31,31,30,31,30,31]
array_ones = np.ones(E_avgmonthly.shape)
for i in range(12):
    array_ones[i] = (days_in_month[i]*array_ones[i])
E_avgmonthly = (E_avgmonthly*(-1000)*(array_ones))
E_avgmonthly_to_year = E_avgmonthly.sum('time')
#E_avgmonthly_to_year.plot()
E_avgmonthly.to_netcdf('/home/chandra/data/Kruger/era5/evap_ERA5_monthly_averaged_reanalysis_2007-resampled_monthly.nc')
