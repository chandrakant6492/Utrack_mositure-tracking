import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy
import numpy as np
from tqdm import tqdm
import pandas as pd

### Add your evaporation data here at 0.5 degree resolution. Mind that it should alighns with the coordinates of moisture recycling data 
### Check the data dimension (360,720), and arrange the (latitude) = (90 to -90), (longitude) = (0, 360)

"""
1. Load the evaporation dataset here
"""
Evap_agg = xr.open_mfdataset('/home/chandra/data/Paper4_Self-influencing_feedback/ERA5_monthly/era5_evaporation_monthly_*_p25.nc').__xarray_dataarray_variable__; # 0.5 degree resolution

# Convert to multi-year mean
Evap_agg = Evap_agg.groupby('time.month').mean(dim = 'time')

# Run this part to convert latitude 90:-90 & longitude 0:360 
"""
lon = Evap_agg.longitude
lon = lon.values
Evap_agg.longitude.values = np.where(lon>= 0,lon,lon+360)
Evap_agg = Evap_agg.sortby(Evap_agg.longitude)
Evap_agg = Evap_agg[:,::-1,:]
Evap_agg
"""

"""
2. Load your dataset for backtrcking here
"""
Data = xr.open_dataset('/home/chandra/data/Kruger/Backtracking/africa_LU-2.nc').get('Africa prop. gridcell') 

# Run this part to convert latitude 90:-90 & longitude 0:360 
lon = Data.longitude
lon = lon.values
Data.longitude.values = np.where(lon>= 0,lon,lon+360)
Data = Data.sortby(Data.longitude)
Data
Screening = (Data/Data.values) # Classifies the sink regions to 1 

lat = []
lon = []
Screened_data = []
for y in range(Screening.shape[0]):
    for x in range(Screening.shape[1]):
        if np.isnan(Screening[y,x].values) == True:
            continue
        else:
            lat.append(round(float(Screening[y,x].latitude.values),2)) #
            lon.append(round(float(Screening[y,x].longitude.values),2))
            Screened_data.append(Screening[y,x].values)

# Extraxt the data as dataframe
Screened_data = pd.DataFrame({'Lat': lat,'Lon': lon,'Screened': Screened_data})

####### Only for trial run 
#Screened_data = Screened_data[0:10]
Screened_data = Screened_data.where(Screened_data['Lat'] == 0.25).dropna()[0:10]

import warnings
warnings.filterwarnings("ignore")

"""
3. Function for moisture recycling
"""
def MR_yearly_source(source_lat,source_lon):
    global main, Evap
    #global main
    main = 0
    month_name = ['01','02','03','04','05','06','07','08','09','10','11','12']
    for i in (range(12)):
        MR = xr.open_dataset('/home/chandra/data/Paper4_Self-influencing_feedback/utrack_climatology/utrack_climatology_0.5_'+
                                         str(month_name[i])+'.nc').moisture_flow
        source = MR.sel(sourcelat = source_lat-0.25, sourcelon= source_lon-0.25)
        source = source.where(source != 255)
        source.values = np.exp(source.values*-0.1)
        Evap = (Evap_agg[i]).sel(lat = slice(source_lat+0.25, source_lat-0.20), 
                                     lon = slice(source_lon-0.20, source_lon+0.25)).values
        if ((np.isnan(Evap) == True) | (Evap == 0)):
        #    print('skipped')
            continue
        else:
            main += source.fillna(0)*Evap
            #print('Recy-total',np.nansum(source*Evap),' ,','Evap',Evap.sum())
        MR.close()
        source.close()

total_main = 0

#### Run in parallel by changing dataframe range 
for i in tqdm(range(Screened_data.shape[0])):
    source_lat, source_lon = np.array(Screened_data[['Lat','Lon']])[i]
    MR_yearly_source(source_lat,source_lon)
    if np.isnan(Evap) == True:
        continue
    else:
        total_main += main

fig = plt.figure(figsize=(7, 5))
ax = [plt.subplot(111,projection=ccrs.PlateCarree(), aspect='auto')]
total_main.plot.contourf(levels = 10, ax = ax[0], transform=ccrs.PlateCarree(), vmin = 0, cmap = 'Reds',)
ax[0].coastlines(lw = 1)
ax[0].add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.3)
ax[0].set_xlim(-20,60)
ax[0].set_ylim(-40, 20)
ax[0].scatter(np.array(Screened_data['Lon']), np.array(Screened_data['Lat']), c='black', alpha = 0.7, s = 0.5, marker = '.')
ax[0].set_title('Forwardtracking (mm/year)')

# Change name here for saving file
total_main.to_netcdf('/home/chandra/data/test_run_forwardtracking.nc')
