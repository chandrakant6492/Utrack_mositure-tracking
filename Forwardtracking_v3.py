import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy
import numpy as np
from tqdm import tqdm, tnrange
import pandas as pd
from scipy.io.netcdf import *
from netCDF4 import Dataset
from numpy import *
import sys

"""
1. Load your evaporation dataset here
"""
print('To run the model, your data format should be:')
print('Latitude = [ 90 → -90 ]')
print('Longitude = [ 0 → 360 ]')
print(' ')

lats=np.arange(90,-90,-0.5)
lons=np.arange(0,360,0.5)

Evap_agg = xr.open_mfdataset('/home/chandra/data/Kruger/era5/evap_ERA5_monthly_averaged_reanalysis_resampled_monthly_mm_per_month_0.5degree.nc').e[:,::-1].sel(time = slice('2008','2017'))
# Convert to multi-year mean
Evap_agg = Evap_agg.groupby('time.month').mean(dim = 'time')
Evap_agg = (Evap_agg.where(Evap_agg > 0))

print('Your Evaporation data:')
print('Latitude = [', Evap_agg.lat[0].values, '→', Evap_agg.lat[-1].values, ']')
print('Longitude = [', Evap_agg.lon[0].values, '→', Evap_agg.lon[-1].values, ']')


"""
2. Load your dataset in Dataframe format
Note: Longitude should be from 0 to 360 degrees 
"""
Screened_data = pd.DataFrame({'Lat': [-18.5, -10.0],'Lon': [310,295]})


### Forward tracking

"""
3. Funtion for forwardtracking (monthly) [Annual can be derived using 'np.sum' funtion] 
"""
def get_closest_index(lats,lat):
        import operator
        lat_index, min_value = min(enumerate(abs(lats-lat)), key=operator.itemgetter(1))
        return lat_index
    
def MR_footprint_forward(month, latitude, longitude):
    latidx=get_closest_index(lats,latitude)
    lonidx=get_closest_index(lons,longitude)
    
    MR = xr.open_dataset('/home/chandra/data/Paper4_Self-influencing_feedback/utrack_climatology/utrack_climatology_0.5_'+
                                                         str(month).zfill(2)+'.nc').moisture_flow
    fp = MR[get_closest_index(lats,latitude), get_closest_index(lons,longitude),:,:].values
    fp=fp*-0.1
    fp=e**fp
    fp[fp==1]=0
    forward_fp=fp/np.nansum(fp)
    #Check (remove comment)
    #print('This values should always be 1 ', np.nansum(forward_fp).round(2))
    forward_fp = forward_fp*((Evap_agg[month-1])[get_closest_index(lats,latitude), 
                                                 get_closest_index(lons,longitude)].values)
    return forward_fp

def absolute_forward_tracking(Screened_data):
    Evap_footprint_monthly = np.zeros(shape=(12,360,720))
    Evap_footprint_monthly_final  = np.zeros(shape=(12,360,720))
    for j in tqdm(range(Screened_data.shape[0])):
        latitude, longitude = np.array(Screened_data.loc[j])
        for i in range(12):
            Evap_footprint_monthly[i,:,:] = MR_footprint_forward(i+1,latitude,longitude)
        Evap_footprint_monthly_final = Evap_footprint_monthly_final + Evap_footprint_monthly

    Evap_footprint_monthly_sum = xr.DataArray(Evap_footprint_monthly_final, coords=[Evap_agg.month.values, Evap_agg.lat.values, Evap_agg.lon.values],
                 dims=['month', 'lat', 'lon'], name = 'Evap_footprint', attrs=dict(description="Evaporation Footprint", units="mm/month"))
    return Evap_footprint_monthly_sum

"""
4. Running the backtracking code for the screened dataframe 
"""
Evap_footprint_monthly_sum = absolute_forward_tracking(Screened_data)


"""
5. Plotting the results
"""
fig = plt.figure(figsize=(3.5, 2.5), dpi = 200)
ax = [plt.subplot(111,projection=ccrs.PlateCarree(), aspect='auto')]
Evap_footprint_monthly_sum.sum(axis = 0).plot(ax = ax[0], transform=ccrs.PlateCarree(), vmin = 0, cmap = 'viridis_r',)
ax[0].coastlines(lw = 1)
ax[0].add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.3)
ax[0].set_xlim(-100,10)
ax[0].set_ylim(-60, 30)
ax[0].scatter(np.array(Screened_data['Lon']), np.array(Screened_data['Lat']), c='red', alpha = 0.7, s = 25, marker = '.')
ax[0].set_title('Forwardtracking: mm/year')


"""
6. Save the dataset 
"""
Evap_footprint_monthly_sum.to_netcdf('/home/chandra/Current/Paper-4/test_forwardtracking.nc')
