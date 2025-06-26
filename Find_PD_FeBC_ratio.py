#%%
import xarray as xr
import pandas as pd
import numpy as np

# FIND MEAN PRESENT DAY FRACTION OF Fe:BC ------------------------------------------------------------------------------
# BC EMISSIONS PD CMIP6
ds_PD = xr.open_dataset("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\BC-em-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_2000-2015MEAN_remapcon_regridded.nc") 
ds_PD = ds_PD.drop_vars(["time_bnds", "sector_bnds", "time"])
ds_PD = ds_PD.squeeze(dim='time') 

# Fe EMISSIONS PD (use as OCNFRAC too)
ds_PD_Fe = xr.open_dataset("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\ClaquinMineralsCAM6_SPEWFUELS-Dec2023_OCNFRAC_remapcon_regridded.nc")
ds_PD_Fe['OCNFRAC'] = ds_PD_Fe['OCNFRAC'].squeeze(dim='time')
ds_PD_Fe = ds_PD_Fe.drop_dims(["time"]) 

# Extract Black Carbon (BC) emission values SSP370 
BC_em_anthro_PD = ds_PD['BC_em_anthro'].copy() # Emission data
lon = ds_PD['lon'].values  # Longitude values -- the .values extracts as a numPy array, removes metadata
lat = ds_PD['lat'].values  # Latitude values
sector = ds_PD['sector'].values # Sector identifiers

# Adding cell area to dataframes so that I can normalize emissions to grid cell area
lon_rads = np.radians(lon)
lat_rads = np.radians(lat)
d_lat = np.abs(lat[1] - lat[0])  # Latitude grid spacing in degrees
d_lon = np.abs(lon[1] - lon[0])  # Longitude grid spacing in degrees
g_lat = np.radians(d_lat / 2)  # Latitude half-spacing in radians
g_lon = np.radians(d_lon / 2)  # Longitude half-spacing in radians

R = 6.3781E6 
cell_areas_staggered = []
for i in range(len(lat)):
    for j in range(len(lon)):
        lat_center = lat_rads[i]
        lon_center = lon_rads[j]
        lat_north = lat_center + g_lat
        lat_south = lat_center - g_lat
        lat_north = np.clip(lat_north, -np.pi / 2, np.pi / 2)
        lat_south = np.clip(lat_south, -np.pi / 2, np.pi / 2)
        area = R**2 * (np.sin(lat_north) - np.sin(lat_south)) * (2 * g_lon)
        cell_areas_staggered.append(area)

cell_areas_staggered = np.array(cell_areas_staggered).reshape(len(lat), len(lon))

# Verify to see if areas add to 5.1E14 in m^2
sum_sa_earth = cell_areas_staggered.sum()
print(f"surface area, {sum_sa_earth:3e}") 

# add cell area to original xarray for easy calling 
ds_PD['cell_area'] = xr.DataArray(
    cell_areas_staggered,
    dims=["lat", "lon"],  # Same dimensions as in the original dataset
    coords={"lat": ds_PD['lat'], "lon": ds_PD['lon']},  # Use original coordinates
    attrs={
        "units": "m^2",  # Specify units for the cell area
        "description": "Calculated grid cell area using staggered grid approach",},)

ds_PD_Fe['cell_area'] = xr.DataArray(
    cell_areas_staggered,
    dims=["lat", "lon"],  # Same dimensions as in the original dataset
    coords={"lat": ds_PD_Fe['lat'], "lon": ds_PD_Fe['lon']},  # Use original coordinates
    attrs={
        "units": "m^2",  # Specify units for the cell area
        "description": "Calculated grid cell area using staggered grid approach",},)

BC_em_anthro_PD['cell_area'] = xr.DataArray(
    cell_areas_staggered,
    dims=["lat", "lon"],  # Same dimensions as in the original dataset
    coords={"lat": BC_em_anthro_PD['lat'], "lon": BC_em_anthro_PD['lon']},  # Use original coordinates
    attrs={
        "units": "m^2",  # Specify units for the cell area
        "description": "Calculated grid cell area using staggered grid approach",},)

# Adding OCNFRAC to separate shipping and terrestrial transportation emissions -------------------------------------------
ocnfrac = ds_PD_Fe['OCNFRAC'] 
ocnfrac_expanded_PD = ocnfrac.expand_dims(dim={'sector': BC_em_anthro_PD['sector']}, axis=0)
# Add OCNFRAC to the BC_emiss dataset
BC_em_anthro_PD['OCNFRAC'] = ocnfrac_expanded_PD
BC_emiss_3_PD = BC_em_anthro_PD.isel(sector=3).copy()  # Transportation
filtered_BC_emiss_3_PD = BC_emiss_3_PD.where(ocnfrac < 0.5)
filtered_BC_emiss_3_PD = filtered_BC_emiss_3_PD.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')
BC_emiss_7_PD = BC_em_anthro_PD.isel(sector=7).copy()  # Shipping
filtered_BC_emiss_7_PD = BC_emiss_7_PD.where(ocnfrac >= 0.5)
filtered_BC_emiss_7_PD = filtered_BC_emiss_7_PD.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
filtered_BC_emiss_3_PD.plot(subplot_kws={'projection': ccrs.PlateCarree()},transform=ccrs.PlateCarree(), cbar_kwargs={'label': 'Emissions'})
plt.title("TRANSPORT")
plt.gca().coastlines()
plt.show() 
filtered_BC_emiss_7_PD.plot(subplot_kws={'projection': ccrs.PlateCarree()},transform=ccrs.PlateCarree(), cbar_kwargs={'label': 'Emissions'})
plt.title("SHIPPING")
plt.gca().coastlines()
plt.show() 
# QA/QC - confirmed that ocean and terrestrial cells were split --------------------------------------------------------

# assigning new arrays for each sector without referencing original xarray --------------------------------------------
BC_emiss_1_PD = BC_em_anthro_PD.isel(sector=1).copy() 
BC_emiss_1_PD = BC_emiss_1_PD.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')# Energy
BC_emiss_2_PD = BC_em_anthro_PD.isel(sector=2).copy()  # Industrial Coal
BC_emiss_2_PD = BC_emiss_2_PD.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')# Energy
BC_emiss_4_PD = BC_em_anthro_PD.isel(sector=4).copy()  # Residential Coal
BC_emiss_4_PD = BC_emiss_4_PD.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')# Energy

CoaMed_FCOALBF_PD = ds_PD_Fe['CoaMed_FCOALBF'].copy() 
FineMed_FCOALBF_PD = ds_PD_Fe['FineMed_FCOALBF'].copy() 
CoaMed_FCOALFF_PD = ds_PD_Fe['CoaMed_FCOALFF'].copy() 
FineMed_FCOALFF_PD = ds_PD_Fe['FineMed_FCOALFF'].copy()   
CoaMed_FOIL_PD = ds_PD_Fe['CoaMed_FOIL'].copy()  
FineMed_FOIL_PD = ds_PD_Fe['FineMed_FOIL'].copy()   
CoaMed_FWOOD_PD = ds_PD_Fe['CoaMed_FWOOD'].copy()  
FineMed_FWOOD_PD = ds_PD_Fe['FineMed_FWOOD'].copy()  

cell_area_array = ds_PD_Fe['cell_area'].copy()

# Finding mass emissions for each grid cell normalized to area -------------------------------------------------------------------
BC_emiss_1_PD_mass = (BC_emiss_1_PD * cell_area_array)/sum_sa_earth
BC_emiss_2_PD_mass = (BC_emiss_2_PD * cell_area_array)/sum_sa_earth
filtered_BC_emiss_3_PD_mass = (filtered_BC_emiss_3_PD * cell_area_array)/sum_sa_earth
BC_emiss_4_PD_mass = (BC_emiss_4_PD * cell_area_array)/sum_sa_earth
filtered_BC_emiss_7_PD_mass = (filtered_BC_emiss_7_PD * cell_area_array)/sum_sa_earth

CoaMed_FCOALBF_PD_mass = (CoaMed_FCOALBF_PD * cell_area_array)/sum_sa_earth
FineMed_FCOALBF_PD_mass = (FineMed_FCOALBF_PD * cell_area_array)/sum_sa_earth
CoaMed_FCOALFF_PD_mass = (CoaMed_FCOALFF_PD * cell_area_array)/sum_sa_earth
FineMed_FCOALFF_PD_mass = (FineMed_FCOALFF_PD * cell_area_array)/sum_sa_earth
CoaMed_FOIL_PD_mass = (CoaMed_FOIL_PD * cell_area_array)/sum_sa_earth  
FineMed_FOIL_PD_mass = (FineMed_FOIL_PD * cell_area_array)/sum_sa_earth
CoaMed_FWOOD_PD_mass = (CoaMed_FWOOD_PD * cell_area_array)/sum_sa_earth
FineMed_FWOOD_PD_mass = (FineMed_FWOOD_PD * cell_area_array)/sum_sa_earth

BC_emiss_1_PD_mass = np.where(BC_emiss_1_PD_mass==0.0, np.nan, BC_emiss_1_PD_mass)
BC_emiss_2_PD_mass = np.where(BC_emiss_2_PD_mass == 0.0, np.nan, BC_emiss_2_PD_mass)
filtered_BC_emiss_3_PD_mass = np.where(filtered_BC_emiss_3_PD_mass == 0.0, np.nan, filtered_BC_emiss_3_PD_mass)
BC_emiss_4_PD_mass = np.where(BC_emiss_4_PD_mass == 0.0, np.nan, BC_emiss_4_PD_mass)
filtered_BC_emiss_7_PD_mass = np.where(filtered_BC_emiss_7_PD_mass == 0.0, np.nan, filtered_BC_emiss_7_PD_mass)

CoaMed_FCOALBF_PD_mass = np.where(CoaMed_FCOALBF_PD_mass == 0.0, np.nan, CoaMed_FCOALBF_PD_mass)
FineMed_FCOALBF_PD_mass = np.where(FineMed_FCOALBF_PD_mass == 0.0, np.nan, FineMed_FCOALBF_PD_mass)
CoaMed_FCOALFF_PD_mass = np.where(CoaMed_FCOALFF_PD_mass == 0.0, np.nan, CoaMed_FCOALFF_PD_mass)
FineMed_FCOALFF_PD_mass = np.where(FineMed_FCOALFF_PD_mass == 0.0, np.nan, FineMed_FCOALFF_PD_mass)
CoaMed_FOIL_PD_mass = np.where(CoaMed_FOIL_PD_mass == 0.0, np.nan, CoaMed_FOIL_PD_mass)
FineMed_FOIL_PD_mass = np.where(FineMed_FOIL_PD_mass == 0.0, np.nan, FineMed_FOIL_PD_mass)
CoaMed_FWOOD_PD_mass = np.where(CoaMed_FWOOD_PD_mass == 0.0, np.nan, CoaMed_FWOOD_PD_mass)
FineMed_FWOOD_PD_mass = np.where(FineMed_FWOOD_PD_mass == 0.0, np.nan, FineMed_FWOOD_PD_mass)

# Total all emissions to calculate global mean
combined_set = np.array([BC_emiss_1_PD_mass, BC_emiss_2_PD_mass, BC_emiss_4_PD_mass, filtered_BC_emiss_3_PD_mass, filtered_BC_emiss_7_PD_mass])
Total_BC_emiss_avg = np.nanmean(combined_set)
combined_set_2 = np.array([CoaMed_FCOALBF_PD_mass, FineMed_FCOALBF_PD_mass, CoaMed_FCOALFF_PD_mass, FineMed_FCOALFF_PD_mass, CoaMed_FOIL_PD_mass, FineMed_FOIL_PD_mass, CoaMed_FWOOD_PD_mass, FineMed_FWOOD_PD_mass])
Total_Fe_emiss_avg = np.nanmean(combined_set_2)

print("BC", Total_BC_emiss_avg)
print("Fe", Total_Fe_emiss_avg)

mean_percent_Fe = Total_Fe_emiss_avg/Total_BC_emiss_avg
print(mean_percent_Fe , "Avg % Iron")

print(CoaMed_FOIL_PD_mass)

# this isnt working right now because OCNFRAC was removed and needs to be indexed instead 
ocnfrac = ds_PD_Fe['OCNFRAC'].values  # Convert to NumPy array if it's not already
terrestrial_mask = ocnfrac < 0.5
oceanic_mask = ocnfrac >= 0.5
terr_CoaMed_FOIL_PD_mass = np.where(terrestrial_mask, CoaMed_FOIL_PD_mass, np.nan)
terr_FineMed_FOIL_PD_mass = np.where(terrestrial_mask, FineMed_FOIL_PD_mass, np.nan)

ocn_CoaMed_FOIL_PD_mass = np.where(oceanic_mask, CoaMed_FOIL_PD_mass, np.nan)
ocn_FineMed_FOIL_PD_mass = np.where(oceanic_mask, FineMed_FOIL_PD_mass, np.nan)

combined_set_3 = np.array([BC_emiss_1_PD_mass, BC_emiss_2_PD_mass, BC_emiss_4_PD_mass, filtered_BC_emiss_3_PD_mass])
terr_BC_emiss_avg = np.nanmean(combined_set_3)
combined_set_4 = np.array([filtered_BC_emiss_7_PD_mass])
ocn_BC_emiss_avg = np.nanmean(combined_set_4)

combined_set_5 = np.array([CoaMed_FCOALBF_PD_mass, FineMed_FCOALBF_PD_mass, CoaMed_FCOALFF_PD_mass, FineMed_FCOALFF_PD_mass, CoaMed_FWOOD_PD_mass, FineMed_FWOOD_PD_mass, terr_CoaMed_FOIL_PD_mass, terr_FineMed_FOIL_PD_mass])
terr_Fe_emiss_avg = np.nanmean(combined_set_5)
combined_set_6 = np.array([ocn_FineMed_FOIL_PD_mass, ocn_CoaMed_FOIL_PD_mass])
ocn_Fe_emiss_avg = np.nanmean(combined_set_6)

print("terr", terr_Fe_emiss_avg)
print("ocn", ocn_Fe_emiss_avg)

mean_percent_Fe_terr = terr_Fe_emiss_avg/terr_BC_emiss_avg
mean_percent_Fe_ocn = ocn_Fe_emiss_avg/ocn_BC_emiss_avg

print(mean_percent_Fe, "Avg % Iron")
print(mean_percent_Fe_terr, "Avg % Iron Terrestrial")
print(mean_percent_Fe_ocn, "Avg % Iron Ocean")

#%%

# Saving mean Fe Fraction for Terrestrial and Marine cells -------------------------------------------------
terr_ratio_Fe_BC_PD = terr_Fe_emiss_sum/terr_BC_emiss_sum
print(terr_ratio_Fe_BC_PD, "Terrestrial % Fe")

ocn_ratio_Fe_BC_PD = ocn_Fe_emiss_sum/ocn_BC_emiss_sum
print(ocn_ratio_Fe_BC_PD, "Ocean % Fe")

# Calculate the average for each
avg_BC_emiss_1_PD_mass = np.nanmean(BC_emiss_1_PD_mass)
avg_BC_emiss_2_PD_mass = np.nanmean(BC_emiss_2_PD_mass)
avg_filtered_BC_emiss_3_PD_mass = np.nanmean(filtered_BC_emiss_3_PD_mass)
avg_BC_emiss_4_PD_mass = np.nanmean(BC_emiss_4_PD_mass)
avg_filtered_BC_emiss_7_PD_mass = np.nanmean(filtered_BC_emiss_7_PD_mass)

avg_CoaMed_FCOALBF_PD_mass = np.nanmean(CoaMed_FCOALBF_PD_mass)
avg_FineMed_FCOALBF_PD_mass = np.nanmean(FineMed_FCOALBF_PD_mass)
avg_CoaMed_FCOALFF_PD_mass = np.nanmean(CoaMed_FCOALFF_PD_mass)
avg_FineMed_FCOALFF_PD_mass = np.nanmean(FineMed_FCOALFF_PD_mass)
avg_CoaMed_FOIL_PD_mass = np.nanmean(CoaMed_FOIL_PD_mass)
avg_FineMed_FOIL_PD_mass = np.nanmean(FineMed_FOIL_PD_mass)
avg_CoaMed_FWOOD_PD_mass = np.nanmean(CoaMed_FWOOD_PD_mass)
avg_FineMed_FWOOD_PD_mass = np.nanmean(FineMed_FWOOD_PD_mass)

# generating sum of BC and Fe emissions from all global cells for each sector -----------------------------------------------------
sum_BC_emiss_1_PD = BC_emiss_1_PD.sum().item()
sum_BC_emiss_2_PD = BC_emiss_2_PD.sum().item()
sum_BC_emiss_3_PD = BC_emiss_4_PD.sum().item()
sum_BC_emiss_4_PD = filtered_BC_emiss_3_PD.sum().item()
sum_BC_emiss_7_PD = filtered_BC_emiss_7_PD.sum().item()

sum_CoaMed_FCOALBF_PD = CoaMed_FCOALBF_PD.sum().item()
sum_FineMed_FCOALBF_PD = FineMed_FCOALBF_PD.sum().item()
sum_CoaMed_FCOALFF_PD = CoaMed_FCOALFF_PD.sum().item()
sum_FineMed_FCOALFF_PD = FineMed_FCOALFF_PD.sum().item()
sum_CoaMed_FOIL_PD = CoaMed_FOIL_PD.sum().item()
sum_FineMed_FOIL_PD = FineMed_FOIL_PD.sum().item()
sum_CoaMed_FWOOD_PD = CoaMed_FWOOD_PD.sum().item()
sum_FineMed_FWOOD_PD = FineMed_FWOOD_PD.sum().item()

# Summing all sectors and calculating ratio for entire globe  ----------------------------------------------------
Total_BC_emiss_sum = (sum_BC_emiss_1_PD + sum_BC_emiss_2_PD + sum_BC_emiss_3_PD + sum_BC_emiss_4_PD + sum_BC_emiss_7_PD)
Total_Fe_emiss_sum = (sum_CoaMed_FCOALBF_PD + sum_FineMed_FCOALBF_PD + sum_CoaMed_FCOALFF_PD + sum_FineMed_FCOALFF_PD + sum_CoaMed_FOIL_PD + sum_FineMed_FOIL_PD + sum_CoaMed_FWOOD_PD + sum_FineMed_FWOOD_PD)
ratio_Fe_BC_PD = (Total_Fe_emiss_sum/Total_BC_emiss_sum)
print(ratio_Fe_BC_PD, "Total % Fe")

# Separating marine emissions from terrestrial emissions  ----------------------------------------------------
terr_CoaMed_FOIL_PD = CoaMed_FOIL_PD.where(ocnfrac < 0.5)
terr_CoaMed_FOIL_PD = terr_CoaMed_FOIL_PD.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')
sum_terr_CoaMed_FOIL_PD = terr_CoaMed_FOIL_PD.sum().item()
terr_FineMed_FOIL_PD = FineMed_FOIL_PD.where(ocnfrac < 0.5)
terr_FineMed_FOIL_PD = terr_CoaMed_FOIL_PD.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')
sum_terr_FineMed_FOIL_PD = terr_FineMed_FOIL_PD.sum().item()

ocn_CoaMed_FOIL_PD = CoaMed_FOIL_PD.where(ocnfrac >= 0.5)
ocn_CoaMed_FOIL_PD = ocn_CoaMed_FOIL_PD.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')
sum_ocn_CoaMed_FOIL_PD = ocn_CoaMed_FOIL_PD.sum().item()
ocn_FineMed_FOIL_PD = FineMed_FOIL_PD.where(ocnfrac >= 0.5)
ocn_FineMed_FOIL_PD = ocn_CoaMed_FOIL_PD.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')
sum_ocn_FineMed_FOIL_PD = ocn_FineMed_FOIL_PD.sum().item()

terr_BC_emiss_sum = (sum_BC_emiss_1_PD + sum_BC_emiss_2_PD + sum_BC_emiss_3_PD + sum_BC_emiss_4_PD)
ocn_BC_emiss_sum = (sum_BC_emiss_7_PD)

terr_Fe_emiss_sum = (sum_CoaMed_FCOALBF_PD + sum_FineMed_FCOALBF_PD + sum_CoaMed_FCOALFF_PD + sum_FineMed_FCOALFF_PD + sum_CoaMed_FWOOD_PD + sum_FineMed_FWOOD_PD + sum_ocn_CoaMed_FOIL_PD + sum_ocn_FineMed_FOIL_PD)
ocn_Fe_emiss_sum = (sum_ocn_FineMed_FOIL_PD + sum_ocn_CoaMed_FOIL_PD)

# Saving mean Fe Fraction for Terrestrial and Marine cells -------------------------------------------------
terr_ratio_Fe_BC_PD = terr_Fe_emiss_sum/terr_BC_emiss_sum
print(terr_ratio_Fe_BC_PD, "Terrestrial % Fe")

ocn_ratio_Fe_BC_PD = ocn_Fe_emiss_sum/ocn_BC_emiss_sum
print(ocn_ratio_Fe_BC_PD, "Ocean % Fe")

# ------------- CODE TO SHOW INTERACTIVE ARRAY WITH iNDIVIDUAL CELL ENTRIES --------------------------------
#Total_BC_emiss_sum_df = pd.DataFrame(BC_emiss_4_PD)
#from itables import show
#show(Total_BC_emiss_sum_df, scrollY=True, scrollX=True, maxRows=100, maxColumns=200) 
#print(type(BC_emiss_1_PD), type(BC_emiss_2_PD), type(BC_emiss_4_PD), type(filtered_BC_emiss_3_PD), type(filtered_BC_emiss_7_PD))


