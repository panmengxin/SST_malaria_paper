import numpy as np 
import matplotlib.pyplot as plt
import datetime
import pandas as pd
import math
import xarray as xr
from multiprocessing import Pool
from scipy import stats

location_save='path of the datasets'


# ## 1. correlation_map_function
# Generates correlation maps between SSTA fields and disease incidence (e.g., malaria) anomaly time series.
# 
# **the dataset we provide:**
# - SSTA field data: the sea surface temperature anomaly fields in 1st epiweek of 1999 to 17th epiweek of 2023
    ## For each grid, the SST anomaly (SSTA) is calculated by subtracting the seasonally varying climatology from the raw SST value to remove the seasonal cycle. 
    ##To match the period of the malaria dataset, we reorganized the daily SSTA field into each epidemiological week. 
    ##For each grid, the linear trend of the SSTA in each calendar epidemiological week is removed. 
# - Disease incidence anomaly time series: malaria anomaly time series over entire loreto region during 1st epiweek of 2000 to 14th epiweek of 2023 (the up-to-date malaria data we have)
    ##the malaria anomaly is used, as we focus on the year-to-year variation. It is calculated by subtracting the malaria annual cycle 
    ##(i.e., mean malaria incidence in each calendar epidemiological week during 2000-2022) from the total malaria incidence.

# Note: To improve the correlation robustness, we concatenate the data for every four epidemiological weeks in the correlation analysis, so you can see the "add" dimension in both SST and Malaria data. It concatenates the data in four sequencial epiweeks. 

# **Output:**
# - Correlation map array containing r-values, p-values, slopes, and intercepts
# 
# with correlation_map_function, you can generate correlation maps for SSTA and disease incidences in different calendar weeks and with different time lags.
# for example, in our study 2756 correlation maps are generated. we group them together and form the "correlation_map_total", as the input of the second function "SOM_function" function.
# Multiple correlation maps can be combined into a `correlation_map_total` array for subsequent clustering analysis.
# correlation maps with longitude and latitude need to be flattened into one dimension (the grids dimension in SOM_function)


# 4 epi week concatenate lead lag correlation map 

SST_anomaly_yearepiweek=xr.open_dataarray(location_save+"SST_anomaly_detrended_40N40S_199901to202317_4epiweekconcatenate.nc")
Malaria_anomaly=xr.open_dataarray(location_save+"Peru_Malaria_anomaly_2000epiweek1_to_2023epiweek14_without53_YearEpiweek_4epiweekconcatenate.nc")
variable_lat=SST_anomaly_yearepiweek['latitude']
variable_long=SST_anomaly_yearepiweek['longitude']


def correlation_map(lead_lag_week):
    correlation_map_xr=xr.DataArray(np.nan,dims=["correlation","epiweek","latitude","longitude"], \
                                    coords={"correlation":["r_value","p_value","slope","intercept"], \
                                            "epiweek":np.arange(1,53),"latitude":variable_lat,"longitude":variable_long})

    # for malaria_epiweek in np.arange(51,52):
    for malaria_epiweek in np.arange(1,53):
        SST_epiweek_old=malaria_epiweek-lead_lag_week
        if malaria_epiweek <= 11:
            if SST_epiweek_old>0:
                SST_epiweek=SST_epiweek_old
                Malaria_anomaly_timeseries=Malaria_anomaly.loc[dict(epiweek=malaria_epiweek,year=slice(2000,2023))]
                SST_anomaly_field=SST_anomaly_yearepiweek.loc[dict(epiweek=SST_epiweek,year=slice(2000,2023))]
            else:
                SST_epiweek=SST_epiweek_old+52
                Malaria_anomaly_timeseries=Malaria_anomaly.loc[dict(epiweek=malaria_epiweek,year=slice(2000,2023))]
                SST_anomaly_field=SST_anomaly_yearepiweek.loc[dict(epiweek=SST_epiweek,year=slice(1999,2022))]
        else:
            if SST_epiweek_old>0:
                SST_epiweek=SST_epiweek_old
                Malaria_anomaly_timeseries=Malaria_anomaly.loc[dict(epiweek=malaria_epiweek,year=slice(2000,2022))]
                SST_anomaly_field=SST_anomaly_yearepiweek.loc[dict(epiweek=SST_epiweek,year=slice(2000,2022))]
            else:
                SST_epiweek=SST_epiweek_old+52
                Malaria_anomaly_timeseries=Malaria_anomaly.loc[dict(epiweek=malaria_epiweek,year=slice(2000,2022))]
                SST_anomaly_field=SST_anomaly_yearepiweek.loc[dict(epiweek=SST_epiweek,year=slice(1999,2021))]            

        print('Malaria anomaly: '+str(malaria_epiweek)+" SSTA: "+str(SST_epiweek))

        add=0
        Malaria_anomaly_timeseries_concatenate=Malaria_anomaly_timeseries.loc[dict(add=add)]
        Malaria_anomaly_timeseries_concatenate=Malaria_anomaly_timeseries_concatenate.assign_coords(year=Malaria_anomaly_timeseries_concatenate['year'].values*10+add)
        SST_anomaly_field_concatenate=SST_anomaly_field.loc[dict(add=add)]
        SST_anomaly_field_concatenate=SST_anomaly_field_concatenate.assign_coords(year=SST_anomaly_field_concatenate['year'].values*10+add)

        for add in range(1,4):
            Malaria_anomaly_timeseries_concatenate_1=Malaria_anomaly_timeseries.loc[dict(add=add)]
            Malaria_anomaly_timeseries_concatenate_1=Malaria_anomaly_timeseries_concatenate_1.assign_coords(year=Malaria_anomaly_timeseries_concatenate_1['year'].values*10+add)
            SST_anomaly_field_concatenate_1=SST_anomaly_field.loc[dict(add=add)]
            SST_anomaly_field_concatenate_1=SST_anomaly_field_concatenate_1.assign_coords(year=SST_anomaly_field_concatenate_1['year'].values*10+add)            
            SST_anomaly_field_concatenate=xr.concat([SST_anomaly_field_concatenate,SST_anomaly_field_concatenate_1],dim='year')
            Malaria_anomaly_timeseries_concatenate=xr.concat([Malaria_anomaly_timeseries_concatenate,Malaria_anomaly_timeseries_concatenate_1],dim='year')

        for lat_index in range(len(variable_lat)):
            lat=variable_lat[lat_index]
            for long_index in range(len(variable_long)):
                long=variable_long[long_index]
                X_series=SST_anomaly_field_concatenate.loc[dict(latitude=lat,longitude=long)].values
                slope, intercept, r_value, p_value, std_err = stats.linregress(X_series,Malaria_anomaly_timeseries_concatenate)
                correlation_map_xr.loc[dict(epiweek=malaria_epiweek,correlation='r_value',latitude=lat,longitude=long)]=r_value
                correlation_map_xr.loc[dict(epiweek=malaria_epiweek,correlation='p_value',latitude=lat,longitude=long)]=p_value
                correlation_map_xr.loc[dict(epiweek=malaria_epiweek,correlation='slope',latitude=lat,longitude=long)]=slope
                correlation_map_xr.loc[dict(epiweek=malaria_epiweek,correlation='intercept',latitude=lat,longitude=long)]=intercept

    correlation_map_xr.to_netcdf(location_save+"correlation_map_SSTA_vs_malaria_anomaly_eachepiweek_lead"+str(lead_lag_week)+"week.nc")
    return lead_lag_week

pool = Pool(20)
parallel=pool.map(correlation_map,np.arange(0,53)) # this is one year


for lead_lag_week in np.arange(0,53):
    correlation_map_xr=xr.open_dataarray(location_save+"correlation_map_SSTA_vs_malaria_anomaly_eachepiweek_lead"+str(lead_lag_week)+"week.nc")
    correlation_map_xr=correlation_map_xr.assign_coords(lead_lag_week=lead_lag_week)
    if lead_lag_week==0:
        correlation_map_xr_total=correlation_map_xr
    else:
        correlation_map_xr_total=xr.concat([correlation_map_xr_total,correlation_map_xr],dim="lead_lag_week")
correlation_map_xr_total.to_netcdf(location_save+"correlation_map_SSTA_vs_malaria_anomaly_eachepiweek_lead_0to52_week.nc")



correlation_map_xr_total=xr.open_dataarray(location_save+"correlation_map_SSTA_vs_malaria_anomaly_eachepiweek_lead_0to52_week.nc").loc[dict(correlation='r_value')]

correlation_map_xr_total_group_long=correlation_map_xr_total.groupby_bins('longitude',np.arange(0,362,2),right=False).mean()
correlation_map_xr_total_group_long=correlation_map_xr_total_group_long.rename(longitude_bins='longitude')
correlation_map_xr_total_group_long=correlation_map_xr_total_group_long.assign_coords(longitude=np.arange(0,360,2))

correlation_map_xr_total_group_long_lat=correlation_map_xr_total_group_long.groupby_bins('latitude',np.arange(-40,42,2),right=False).mean()
correlation_map_xr_total_group_long_lat=correlation_map_xr_total_group_long_lat.rename(latitude_bins='latitude')
correlation_map_xr_total_group_long_lat=correlation_map_xr_total_group_long_lat.assign_coords(latitude=np.arange(-40,40,2))
correlation_map_xr_total_group_long_lat.to_netcdf(location_save+"correlation_map_SSTA_vs_malaria_anomaly_eachepiweek_2degree_lead_0to52_week.nc")

correlation_map_xr_total=xr.open_dataarray(location_save+"Malaria/correlation_map_SSTA_vs_malaria_anomaly_eachepiweek_2degree_lead_0to52_week.nc").loc[dict(correlation='r_value')]

lead_lag_week_list=correlation_map_xr_total['lead_lag_week'].values
epiweek_list=correlation_map_xr_total['epiweek'].values
variable_lat=correlation_map_xr_total['latitude'].values
variable_long=correlation_map_xr_total['longitude'].values

epiweek_leadtime_list=[]
for lead_lag_week in lead_lag_week_list:
    for epiweek in epiweek_list:
        epiweek_leadtime=epiweek*1000+lead_lag_week
        epiweek_leadtime_list.append(epiweek_leadtime)
        
reorganized_correlation_map_xr=xr.DataArray(np.nan,dims=["epiweek_leadtime","latitude","longitude"], \
                                            coords={"epiweek_leadtime":epiweek_leadtime_list,"latitude":variable_lat,"longitude":variable_long})

epiweek_leadtime_list=[]
for lead_lag_week in lead_lag_week_list:
    for epiweek in epiweek_list:
        epiweek_leadtime=epiweek*1000+lead_lag_week
        reorganized_correlation_map_xr.loc[dict(epiweek_leadtime=epiweek_leadtime)]=correlation_map_xr_total.loc[dict(epiweek=epiweek,lead_lag_week=lead_lag_week)]

reorganized_correlation_map_xr.to_netcdf(location_save+"correlation_map_SSTA_vs_malaria_anomaly_eachepiweek_2degree_lead_0to52_week_reorganized_forSOM.nc")


reorganized_correlation_map_xr=xr.open_dataarray(location_save+"correlation_map_SSTA_vs_malaria_anomaly_eachepiweek_2degree_lead_0to52_week_reorganized_forSOM.nc")
reorganized_correlation_map_xr_large=reorganized_correlation_map_xr.where(abs(reorganized_correlation_map_xr)>0.2,0)

reorganized_correlation_map_xr_example=reorganized_correlation_map_xr[0,:,:]
target_index=np.where(np.isnan(reorganized_correlation_map_xr_example)==False)
lat_index=target_index[0]
lon_index=target_index[1]
num_grid=len(lon_index)
AR_freq_total_01=np.zeros(np.shape(reorganized_correlation_map_xr_example))
AR_freq_total_01[target_index]=1

epiweek_leadtime_list=reorganized_correlation_map_xr_large["epiweek_leadtime"].values
reorganized_correlation_series=reorganized_correlation_map_xr_large.values[:,lat_index,lon_index]
reorganized_correlation_series_forSOM=xr.DataArray(reorganized_correlation_series,dims=["epiweek_leadtime","grid_index"], \
                          coords={"epiweek_leadtime":epiweek_leadtime_list,"grid_index":np.arange(0,num_grid)})

reorganized_correlation_series_forSOM.to_netcdf(location_save+"correlation_map_SSTA_vs_malaria_anomaly_eachepiweek_2degree_lead_0to52_week_reorganized_forSOM_Rmorethan02.nc")

