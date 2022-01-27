import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy
import numpy as np
from tqdm import tqdm
import pandas as pd

### Add your precipitation data here at 0.5 degree resolution. Mind that it should alighns with the coordinates of moisture recycling data 
### Check the data dimension (360,720), and arrange the (latitude) = (90 to -90), (longitude) = (0, 360)

"""
1. Load the precipitation dataset here
"""
Precip_agg = xr.open_mfdataset('/home/chandra/data/Paper4_Self-influencing_feedback/CHIRPS/chirps-v2.0.monthly_*_p50.nc').__xarray_dataarray_variable__.sel(time = slice('2006','2015'))

# Convert to multi-year mean
Precip_agg = Precip_agg.groupby('time.month').mean(dim='time')

# Run this part to convert latitude 90:-90 & longitude 0:360 
lon = Precip_agg.longitude
lon = lon.values
Precip_agg.longitude.values = np.where(lon>= 0,lon,lon+360)
Precip_agg = Precip_agg.sortby(Precip_agg.longitude)
Precip_agg = Precip_agg[:,::-1,:]
Precip_agg

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
            lat.append(round(float(Screening[y,x].latitude.values),2))
            lon.append(round(float(Screening[y,x].longitude.values),2))
            Screened_data.append(Screening[y,x].values)

# Extraxt the data as dataframe
Screened_data = pd.DataFrame({'Lat': lat,'Lon': lon,'Screened': Screened_data})

####### Only for trial run 
Screened_data = Screened_data.where(Screened_data['Lat'] == 0.25).dropna()[0:10]

import warnings
warnings.filterwarnings("ignore")

"""
3. Function for moisture recycling
"""
def MR_yearly_sink(sink_lat,sink_lon):
    global main, Precip
    #global main
    main = 0
    month_name = ['01','02','03','04','05','06','07','08','09','10','11','12']
    for i in (range(12)):
        MR = xr.open_dataset('/home/chandra/data/Paper4_Self-influencing_feedback/utrack_climatology/utrack_climatology_0.5_'+
                                         str(month_name[i])+'.nc').moisture_flow
        target = MR.sel(targetlat = sink_lat-0.25, targetlon= sink_lon-0.25)
        target = target.where(target != 255)
        target.values = np.exp(target.values*-0.1)
        Precip = (Precip_agg[i]).sel(latitude = slice(sink_lat+0.05, sink_lat-0.05), 
                                     longitude = slice(sink_lon-0.05, sink_lon+0.05)).values
        if ((np.isnan(Precip) == True) | (Precip == 0)):
        #    print('skipped')
            continue
        else:
            main += target.fillna(0)*Precip
            #print('Recy-total',np.nansum(target*Precip),' ,','Precip',Precip.sum())
        MR.close()
        target.close()

total_main = 0

#### Run in parallel by changing dataframe range 
for i in tqdm(range(Screened_data.shape[0])):
    sink_lat, sink_lon = np.array(Screened_data[['Lat','Lon']])[i]
    MR_yearly_sink(sink_lat,sink_lon)
    if np.isnan(Precip) == True:
        continue
    else:
        total_main += main

fig = plt.figure(figsize=(7, 5), dpi = 200)
ax = [plt.subplot(111,projection=ccrs.PlateCarree(), aspect='auto')]
total_main.plot.contourf(levels = 10, ax = ax[0], transform=ccrs.PlateCarree(), vmin = 0, cmap = 'Blues',)
ax[0].coastlines(lw = 1)
ax[0].add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.3)
ax[0].set_xlim(-20,60)
ax[0].set_ylim(-40, 20)
ax[0].scatter(np.array(Screened_data['Lon']), np.array(Screened_data['Lat']), c='black', alpha = 0.7, s = 0.5, marker = '.')
ax[0].set_title('mm/year')

# Change name here for saving file
total_main.to_netcdf('/home/chandra/data/Kruger/Backtracking/test_run_backtracking.nc')