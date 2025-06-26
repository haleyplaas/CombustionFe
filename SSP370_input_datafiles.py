#%% CHECKING TO SEE IF REGRIDDER WORKED PROPERLY USING CDO REMAPCON --------------------------------------------------
# PASSED FOR PD FILES -- budgets add to the same value
import numpy as np
import xarray as xr

#OG_PD = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\BC-em-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_2000-2015MEAN_remapcon_regridded.nc"
OG_PD = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\BC-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-AIM-ssp370-1-1_gn_2090_2100MEAN_remapcon_regridded.nc"
ds_OG_PD = xr.open_dataset(OG_PD)

BC_em_anthro_OG = ds_OG_PD['BC_em_anthro']
sector = ds_OG_PD['sector'].values # Sector identifiers
#BC_em_anthro_OG1 = BC_em_anthro_OG.isel(sector=4) 

#sum_OG_PD_1 = BC_em_anthro_OG1.sum()
#print("sum", sum_OG_PD_1)

# actually need to calculate emission budgets to determine if values are the same 
lon = ds_OG_PD['lon'].values  # Longitude in degrees
lat = ds_OG_PD['lat'].values  # Latitude in degrees
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
ds_OG_PD['cell_area'] = xr.DataArray(
    cell_areas_staggered,
    dims=["lat", "lon"],  # Same dimensions as in the original dataset
    coords={"lat": ds_OG_PD['lat'], "lon": ds_OG_PD['lon']},  # Use original coordinates
    attrs={
        "units": "m^2",  # Specify units for the cell area
        "description": "Calculated grid cell area using staggered grid approach",
    },
)

import pandas as pd
# Dictionaries to store both individual and total budgets
global_budgets = {}
individual_cell_emissions = {}

# Loop over all variables in the dataset
for var_name in ds_OG_PD.data_vars:
    if "BC_em_anthro" in var_name:  
        # selecting specific sector 
        sector_data = ds_OG_PD[var_name].sel(sector=7)

        # Calculate individual emissions (for each cell) and convert from sec to annual and from kg to Tg
        individual_cell_emissions[var_name] = (sector_data * ds_OG_PD['cell_area'] * 3600 * 24 * 365 *1E-9)
        
        # Calculate the total budget by summing individual emissions
        total_budget = individual_cell_emissions[var_name].sum().values  # Extract scalar value
        
        # Store the total budget in the dictionary
        global_budgets[var_name] = total_budget
        
        # Print the budget for the variable
        print(f"FU Total Fe budget of {var_name} (Tg/year) for sector: {total_budget:.3e}")

# test to see if BC emissions are supposed to be higher for FU
# 0: Agriculture; 1: Energy;  2: Industrial; 3: Transportation; 4: Residential, Commercial, Other; 5: Solvents production and application; 6: Waste; 7: International Shipping
# FOIL = 3 (only from terrestrial cells) + 7 (only from ocean cells)
# FWOOD = 4 
# FSMELT = bring over from Fe_Emissions_Fuel
# FCOALFF = 1 + 2
# FCOALBF = 4

# emission budget for PD sectors 1:9.814e-01, 2:7.405e-01, 3:1.237e+00, 4:3.448e+00, 7:1.705e-01
# emission budget for FU sectors 1:1.618e+00, 2:4.168e-01, 3:1.528e+00, 4:2.881e+00, 7:4.812e-02
# yes, confirmed that OIL and FCOALFF should see an increase in FU scenario per increase to BC, at least one of two sectors increased 


#%% CALCULATING SECTOR SPECIFIC Fe EMISSIONS FROM CMIP BLACK CARBON EMISSION DATASETS  ----------------------
import xarray as xr
import pandas as pd
import numpy as np

# Load the NetCDF file -- BC EMISSIONS SSP370
# ds_SSP370 = xr.open_dataset("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\BC-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-AIM-ssp370-1-1_gn_2090_2100MEAN_remapcon_regridded.nc")
ds_SSP370 = xr.open_dataset("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\BC-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-AIM-ssp370-1-1_gn_2040-2050MEAN_remapcon_regridded.nc")
ds_SSP370 = ds_SSP370.drop_vars(["time_bnds", "sector_bnds", "time"])
ds_SSP370 = ds_SSP370.squeeze(dim='time') 
# # BOUNDING HAS BEEN AN ISSUE -- despite remapping, final value was 90.00000058 instead of 90.

# Load the NetCDF file -- BC EMISSIONS PD CMIP6
ds_PD = xr.open_dataset("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\BC-em-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_2000-2015MEAN_remapcon_regridded.nc") 
ds_PD = ds_PD.drop_vars(["time_bnds", "sector_bnds", "time"])
#ds_PD = ds_PD.drop_dims(["time"])
ds_PD = ds_PD.squeeze(dim='time') 

# Load the NetCDF file -- Fe EMISSIONS PD CMIP6 (use as OCNFRAC too)
ds_PD_Fe = xr.open_dataset("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\ClaquinMineralsCAM6_SPEWFUELS-Dec2023_OCNFRAC_remapcon_regridded.nc")
ds_PD_Fe['OCNFRAC'] = ds_PD_Fe['OCNFRAC'].squeeze(dim='time')
ds_PD_Fe = ds_PD_Fe.drop_dims(["time"]) 
#ds_PD_Fe = ds_PD_Fe.squeeze(dim='time')

#print("BC Future", ds_SSP370.coords)
#print("BC PD", ds_PD.coords)
#print("Fe PD", ds_PD_Fe.coords)
# remapcon ensured that now they all align other than sector 

# Extract Black Carbon (BC) emission values SSP370 
BC_em_anthro_SSP370 = ds_SSP370['BC_em_anthro'] # Emission data
lon = ds_SSP370['lon'].values  # Longitude values -- the .values extracts as a numPy array, removes metadata
lat = ds_SSP370['lat'].values  # Latitude values
sector = ds_SSP370['sector'].values # Sector identifiers

# 0: Agriculture; 1: Energy;  2: Industrial; 3: Transportation; 4: Residential, Commercial, Other; 5: Solvents production and application; 6: Waste; 7: International Shipping
# FOIL = 3 (only from terrestrial cells) + 7 (only from ocean cells)
# FWOOD = 4 
# FSMELT = bring over from Fe_Emissions_Fuel
# FCOALFF = 1 + 2
# FCOALBF = 4

# Specifying BC emission sector for each dataset SSP370 -------------------------------------------------------------
BC_emiss_1_SSP370 = BC_em_anthro_SSP370.isel(sector=1) # Energy
BC_emiss_2_SSP370 = BC_em_anthro_SSP370.isel(sector=2) # Industrial Coal
BC_emiss_4_SSP370 = BC_em_anthro_SSP370.isel(sector=4) # Residential Coal

# Adding OCNFRAC to separate shipping and terrestrial transportation emissions -------------------------------------------
ocnfrac = ds_PD_Fe['OCNFRAC'] 
ocnfrac_expanded_SSP370 = ocnfrac.expand_dims(dim={'sector': BC_em_anthro_SSP370['sector']}, axis=0)

# Add OCNFRAC to the BC_emiss dataset
BC_em_anthro_SSP370['OCNFRAC'] = ocnfrac_expanded_SSP370
BC_emiss_3_SSP370 = BC_em_anthro_SSP370.isel(sector=3).copy()  # Transportation
filtered_BC_emiss_3_SSP370 = BC_emiss_3_SSP370.where(ocnfrac < 0.5)
filtered_BC_emiss_3_SSP370 = filtered_BC_emiss_3_SSP370.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')

BC_emiss_7_SSP370 = BC_em_anthro_SSP370.isel(sector=7).copy()  # Shipping
filtered_BC_emiss_7_SSP370 = BC_emiss_7_SSP370.where(ocnfrac >= 0.5)
filtered_BC_emiss_7_SSP370 = filtered_BC_emiss_7_SSP370.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')

# --------------- Quick plots to ensure that emission data was extracted properly ----------------------
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

#filtered_BC_emiss_3_SSP370.plot(subplot_kws={'projection': ccrs.PlateCarree()},transform=ccrs.PlateCarree(), cbar_kwargs={'label': 'Emissions'})
#plt.title("TRANSPORT")
#plt.gca().coastlines()
#plt.show() # YAY! The land-ocean masking worked this way, checked both transport and shipping

#filtered_BC_emiss_7_SSP370.plot(subplot_kws={'projection': ccrs.PlateCarree()},transform=ccrs.PlateCarree(), cbar_kwargs={'label': 'Emissions'})
#plt.title("SHIPPING")
#plt.gca().coastlines()
#plt.show() # YAY! The land-ocean masking worked this way, checked both transport and shipping

# need to find ratio of BF to WOOD before assigning Residential coal (4) ----------------------------------------
ds_PD_Fe = ds_PD_Fe.drop_vars(['time'], errors= 'ignore')

CoaMed_FCOALBF_PD = ds_PD_Fe['CoaMed_FCOALBF'] 
FineMed_FCOALBF_PD = ds_PD_Fe['FineMed_FCOALBF']   
CoaMed_FCOALFF_PD = ds_PD_Fe['CoaMed_FCOALFF']  
FineMed_FCOALFF_PD = ds_PD_Fe['FineMed_FCOALFF']  
CoaMed_FOIL_PD = ds_PD_Fe['CoaMed_FOIL'] 
FineMed_FOIL_PD = ds_PD_Fe['FineMed_FOIL']   
CoaMed_FSMELT_PD = ds_PD_Fe['CoaMed_FSMELT'] 
FineMed_FSMELT_PD = ds_PD_Fe['FineMed_FSMELT'] 
CoaMed_FWOOD_PD = ds_PD_Fe['CoaMed_FWOOD']  
FineMed_FWOOD_PD = ds_PD_Fe['FineMed_FWOOD']  


# Adding fractions of Fe together
BothMode_FCOALBF_PD = (CoaMed_FCOALBF_PD + FineMed_FCOALBF_PD)
BothMode_FCOALFF_PD = (CoaMed_FCOALFF_PD + FineMed_FCOALFF_PD)
BothMode_FOIL_PD = (CoaMed_FOIL_PD + FineMed_FOIL_PD)
BothMode_FWOOD_PD = (CoaMed_FWOOD_PD + FineMed_FWOOD_PD)

ResiCoal_BF_frac = BothMode_FCOALBF_PD/(BothMode_FCOALBF_PD + BothMode_FWOOD_PD)
average_ResiCoal_BF_frac = ResiCoal_BF_frac.where(np.isfinite(ResiCoal_BF_frac)).mean()
print("average_ResiCoal_BF_frac", average_ResiCoal_BF_frac)

ResiCoal_WOOD_frac = BothMode_FWOOD_PD/(BothMode_FCOALBF_PD + BothMode_FWOOD_PD)
average_ResiCoal_WOOD_frac = ResiCoal_WOOD_frac.where(np.isfinite(ResiCoal_WOOD_frac)).mean()
print("average_ResiCoal_WOOD_frac", average_ResiCoal_WOOD_frac)

Total_ResiCoal = ResiCoal_BF_frac + ResiCoal_WOOD_frac
average_ResiCoal_TOT_frac = Total_ResiCoal.where(np.isfinite(Total_ResiCoal )).mean()
print("average_ResiCoal_TOT_frac", average_ResiCoal_TOT_frac) # QA/QC adds to 1

# Splitting ResiCoal into FCOALBF and WOOD by ratios 
# FCOALBF --------------------------------------------------------------------
zero_mask_FCOALBF_PD = (BC_emiss_4_SSP370 != 0) & ((BothMode_FCOALBF_PD == 0))
BC_COALBF_SSP370_masked = BC_emiss_4_SSP370.where(zero_mask_FCOALBF_PD)
masked_indices_FCOALBF_PD = list(zip(*np.where(zero_mask_FCOALBF_PD)))
print(len(masked_indices_FCOALBF_PD))
print("indices", masked_indices_FCOALBF_PD)

zero_mask_FWOOD_PD = (BC_emiss_4_SSP370 != 0) & ((BothMode_FWOOD_PD == 0))
BC_WOOD_SSP370_masked = BC_emiss_4_SSP370.where(zero_mask_FWOOD_PD)
masked_indices_FWOOD_PD = list(zip(*np.where(zero_mask_FWOOD_PD)))
print(len(masked_indices_FWOOD_PD))
print("indices", masked_indices_FWOOD_PD)

BC_COALBF_SSP370 = xr.zeros_like(BC_emiss_4_SSP370)
BC_WOOD_SSP370 = xr.zeros_like(BC_emiss_4_SSP370)

BC_emiss_4_SSP370_copy = BC_emiss_4_SSP370.copy()

for index in np.ndindex(BC_emiss_4_SSP370.shape):
    if index in masked_indices_FCOALBF_PD:
        # For masked indices (where no PD data), multiply the cell by the average breakdown of BF to Wood
        BC_COALBF_SSP370[index] = BC_emiss_4_SSP370_copy[index] * average_ResiCoal_BF_frac
    else:
        # For all other indices, use cell by cell ratios of BF to Wood
        BC_COALBF_SSP370[index] = BC_emiss_4_SSP370_copy[index] * ResiCoal_BF_frac[index]

    if index in masked_indices_FWOOD_PD:
        # For masked indices (where no PD data), multiply the cell by the average breakdown of BF to Wood
        BC_WOOD_SSP370[index] = BC_emiss_4_SSP370_copy[index] * average_ResiCoal_WOOD_frac
    else:
        # For all other indices, use cell by cell ratios of BF to Wood
        BC_WOOD_SSP370[index] = BC_emiss_4_SSP370_copy[index] * ResiCoal_WOOD_frac[index]

print(BC_COALBF_SSP370)
print(BC_WOOD_SSP370)

#BC_bf_sum = BC_COALBF_SSP370.sum().item()
#BC_wood_sum = BC_WOOD_SSP370.sum().item()
#Total_ResiCoal_sum = BC_emiss_4_SSP370.sum().item()
#added_sum = BC_bf_sum + BC_wood_sum

#print(BC_bf_sum)
#print(BC_wood_sum)
#print(Total_ResiCoal_sum)
#print(added_sum) # tests passed, only slight difference between sums of split BC and original emiss_4, likely explained by the cells where the average was supplied to avoid missing data

BC_COALBF_SSP370 = BC_COALBF_SSP370.fillna(0.0)
BC_WOOD_SSP370 = BC_WOOD_SSP370.fillna(0.0)
BC_COALFF_SSP370 = (BC_emiss_1_SSP370 + BC_emiss_2_SSP370)
BC_OIL_SSP370 = filtered_BC_emiss_3_SSP370.fillna(filtered_BC_emiss_7_SSP370) # replacing all NAs for transportation with 
# SMELT WILL COME DIRECTLY WHEN I ADD IN THE PD DATA

#BC_COALBF_SSP370.plot(subplot_kws={'projection': ccrs.PlateCarree()},transform=ccrs.PlateCarree(), cbar_kwargs={'label': 'Emissions'})
#plt.title("BIOFUEL")
#plt.gca().coastlines()
#plt.show()

#BC_WOOD_SSP370.plot(subplot_kws={'projection': ccrs.PlateCarree()},transform=ccrs.PlateCarree(), cbar_kwargs={'label': 'Emissions'})
#plt.title("WOOD")
#plt.gca().coastlines()
#plt.show()

# Check combination of shipping and transportation with quick plot
#BC_OIL_SSP370.plot(subplot_kws={'projection': ccrs.PlateCarree()},transform=ccrs.PlateCarree(), cbar_kwargs={'label': 'Emissions'})
#plt.title("BC_OIL_SSP370")
#plt.gca().coastlines()
#plt.show() # plot looks right for each sector

# -----------------------------------------------------------------------------------------------------------
# Extract Black Carbon (BC) emission values PD 
BC_em_anthro_PD = ds_PD['BC_em_anthro']  # Emission data
ocnfrac_expanded_PD = ocnfrac.expand_dims(dim={'sector': BC_em_anthro_PD['sector']}, axis=0)

# Add OCNFRAC to the BC_emiss dataset
BC_em_anthro_PD['OCNFRAC'] = ocnfrac_expanded_PD
BC_emiss_3_PD = BC_em_anthro_PD.isel(sector=3).copy()  # Transportation
filtered_BC_emiss_3_PD = BC_emiss_3_PD.where(ocnfrac < 0.5)
filtered_BC_emiss_3_PD = filtered_BC_emiss_3_PD.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')

BC_emiss_7_PD = BC_em_anthro_PD.isel(sector=7).copy()  # Shipping
filtered_BC_emiss_7_PD = BC_emiss_7_PD.where(ocnfrac >= 0.5)
filtered_BC_emiss_7_PD = filtered_BC_emiss_7_PD.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')

# Specifying emission sector for each dataset
BC_emiss_1_PD = BC_em_anthro_PD.isel(sector=1) # Energy
BC_emiss_2_PD = BC_em_anthro_PD.isel(sector=2) # Industrial Coal
BC_emiss_4_PD = BC_em_anthro_PD.isel(sector=4) # Residential Coal

# Separating Residential emissions into biofuel and wood 
zero_mask_FCOALBF_PD = (BC_emiss_4_PD != 0) & ((BothMode_FCOALBF_PD == 0))
BC_COALBF_PD_masked = BC_emiss_4_PD.where(zero_mask_FCOALBF_PD)
masked_indices_FCOALBF_PD = list(zip(*np.where(zero_mask_FCOALBF_PD)))
print(len(masked_indices_FCOALBF_PD))
print("indices", masked_indices_FCOALBF_PD)

zero_mask_FWOOD_PD = (BC_emiss_4_PD != 0) & ((BothMode_FWOOD_PD == 0))
BC_WOOD_PD_masked = BC_emiss_4_PD.where(zero_mask_FWOOD_PD)
masked_indices_FWOOD_PD = list(zip(*np.where(zero_mask_FWOOD_PD)))
print(len(masked_indices_FWOOD_PD))
print("indices", masked_indices_FWOOD_PD)

BC_COALBF_PD = xr.zeros_like(BC_emiss_4_PD)
BC_WOOD_PD = xr.zeros_like(BC_emiss_4_PD)
BC_emiss_4_PD_copy = BC_emiss_4_PD.copy()
for index in np.ndindex(BC_emiss_4_PD.shape):
    if index in masked_indices_FCOALBF_PD:
        # For masked indices (where no PD data), multiply the cell by the average breakdown of BF to Wood
        BC_COALBF_PD[index] = BC_emiss_4_PD_copy[index] * average_ResiCoal_BF_frac
    else:
        # For all other indices, use cell by cell ratios of BF to Wood
        BC_COALBF_PD[index] = BC_emiss_4_PD_copy[index] * ResiCoal_BF_frac[index]
    if index in masked_indices_FWOOD_PD:
        # For masked indices (where no PD data), multiply the cell by the average breakdown of BF to Wood
        BC_WOOD_PD[index] = BC_emiss_4_PD_copy[index] * average_ResiCoal_WOOD_frac
    else:
        # For all other indices, use cell by cell ratios of BF to Wood
        BC_WOOD_PD[index] = BC_emiss_4_PD_copy[index] * ResiCoal_WOOD_frac[index]

BC_COALBF_PD = BC_COALBF_PD.fillna(0.0)
BC_WOOD_PD = BC_WOOD_PD.fillna(0.0)
BC_COALFF_PD = (BC_emiss_1_PD + BC_emiss_2_PD)
BC_OIL_PD = filtered_BC_emiss_3_PD.fillna(filtered_BC_emiss_7_PD) 

#BC_bf_sum = BC_COALBF_PD.sum().item()
#BC_wood_sum = BC_WOOD_PD.sum().item()
#Total_ResiCoal_sum = BC_emiss_4_PD.sum().item()
#added_sum = BC_bf_sum + BC_wood_sum 

#print(BC_bf_sum)
#print(BC_wood_sum)
#print(Total_ResiCoal_sum)
#print(added_sum) # tests passed, only slight difference between sums of split BC and original emiss_4, likely explained by the cells where the average was supplied to avoid missing data

# Check combination 
BC_WOOD_PD.plot(subplot_kws={'projection': ccrs.PlateCarree()},transform=ccrs.PlateCarree(), cbar_kwargs={'label': 'Emissions'})
plt.title("BC_OIL_PD")
plt.gca().coastlines()
plt.show() # plot looks right for each sector

BC_COALBF_PD.plot(subplot_kws={'projection': ccrs.PlateCarree()},transform=ccrs.PlateCarree(), cbar_kwargs={'label': 'Emissions'})
plt.title("BC_COALBF_PD")
plt.gca().coastlines()
plt.show() # plot looks right for each sector

#---------------------------------------------------------------------------------------------------------------
# Finding average ratio of Fe to BC normalized to fuel type mass fractions
total_Fe_emiss = (CoaMed_FCOALBF_PD + FineMed_FCOALBF_PD + CoaMed_FCOALFF_PD + FineMed_FCOALFF_PD + CoaMed_FOIL_PD + FineMed_FOIL_PD + CoaMed_FWOOD_PD + FineMed_FWOOD_PD) # + CoaMed_FSMELT_PD + FineMed_FSMELT_PD
# think I need to not include smelting because we dont have those BC emissions included which is inflating the % 

total_BC_emiss = (BC_COALFF_PD + BC_OIL_PD + BC_COALBF_PD + BC_WOOD_PD)

mask = total_BC_emiss != 0

# Apply the mask to each dataset
masked_ds1 = total_Fe_emiss.where(mask)
masked_ds2 = total_BC_emiss.where(mask)

avg_iron_frac = (total_Fe_emiss/total_BC_emiss)
#avg_iron_frac = (masked_ds1/masked_ds2)

avg_iron_frac_valid = np.ma.masked_invalid(avg_iron_frac).mean()
#avg_iron_frac_valid = avg_iron_frac.mean()

print("avg Fe %", avg_iron_frac_valid)
# the same whether I apply the mask or just ignore the NaNs 

# Masking by OCNFRAC and finding % Fe ratios
# ocean oil = filtered_BC_emiss_7_PD
ocean_mask = ~filtered_BC_emiss_7_PD.isnull() # ~ is is NOT in Python (as opposed to !=)

total_Fe_emiss_OIL = CoaMed_FOIL_PD + FineMed_FOIL_PD 
ocean_Fe_emiss_OIL = total_Fe_emiss_OIL.where(ocean_mask)

ocean_Fe_frac = ocean_Fe_emiss_OIL/filtered_BC_emiss_7_PD
avg_iron_frac_ocean = np.ma.masked_invalid(ocean_Fe_frac).mean()
print("ocean avg Fe %", avg_iron_frac_ocean)

# removing variables I no longer need (OCNFRAC and time from emission xarrays)
# These may need to be removed at the very end 
BC_COALFF_SSP370 = BC_COALFF_SSP370.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')
BC_COALBF_SSP370 = BC_COALBF_SSP370.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')
BC_WOOD_SSP370 = BC_WOOD_SSP370.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')
BC_OIL_SSP370 = BC_OIL_SSP370.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')

BC_COALFF_PD = BC_COALFF_PD.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')
BC_COALBF_PD = BC_COALBF_PD.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')
BC_WOOD_PD = BC_WOOD_PD.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')
BC_OIL_PD = BC_OIL_PD.drop_vars(['time','OCNFRAC','sector'], errors= 'ignore')

# Masks for indexing cells where a BC value exists from SSP370 dataset, but no data is available for PD iron, resulting in an extra NA value
zero_mask_COALBF = (BC_COALBF_SSP370 != 0) & ((BothMode_FCOALBF_PD == 0) | (BC_COALBF_PD == 0))
BC_COALBF_SSP370_masked = BC_COALBF_SSP370.where(zero_mask_COALBF)
masked_indices_COALBF = list(zip(*np.where(zero_mask_COALBF)))
print(len(masked_indices_COALBF))
print("indices", masked_indices_COALBF)

zero_mask_COALFF = (BC_COALFF_SSP370 != 0) & ((BothMode_FCOALFF_PD == 0) | (BC_COALFF_PD == 0))
BC_COALFF_SSP370_masked = BC_COALFF_SSP370.where(zero_mask_COALFF)
masked_indices_COALFF = list(zip(*np.where(zero_mask_COALFF)))
print(len(masked_indices_COALFF))
print("indices", masked_indices_COALFF)

zero_mask_OIL = (BC_OIL_SSP370 != 0) & ((BothMode_FOIL_PD == 0) | (BC_OIL_PD == 0))
BC_OIL_SSP370_masked = BC_OIL_SSP370.where(zero_mask_OIL)
masked_indices_OIL = list(zip(*np.where(zero_mask_OIL)))
print(len(masked_indices_OIL))
print("indices", masked_indices_OIL)

zero_mask_WOOD = (BC_WOOD_SSP370 != 0) & ((BothMode_FWOOD_PD == 0) | (BC_WOOD_PD == 0))
BC_WOOD_SSP370_masked = BC_WOOD_SSP370.where(zero_mask_WOOD)
masked_indices_WOOD = list(zip(*np.where(zero_mask_WOOD)))
print(len(masked_indices_WOOD))
print("indices", masked_indices_WOOD)

# Calculating Fe_SSP370 emission values from ratio of FePD/BCPD to BCFU ----------------------------
# FCOALBF --------------------------------------------------------------------
BothMode_FCOALBF_SSP370 = xr.zeros_like(BC_COALBF_SSP370)
BC_COALBF_SSP370_copy = BC_COALBF_SSP370.copy()
for index in np.ndindex(BC_COALBF_SSP370.shape):
    if index in masked_indices_COALBF:
        # For masked indices (where no PD data), multiply the cell by avg_iron_frac_valid 
        BothMode_FCOALBF_SSP370[index] = BC_COALBF_SSP370_copy[index] * avg_iron_frac_valid
    else:
        # For all other indices, use cell by cell ratios of Fe:BC to determine regional Fe emission frac 
        BothMode_FCOALBF_SSP370[index] = (BothMode_FCOALBF_PD[index] / BC_COALBF_PD[index]) * BC_COALBF_SSP370_copy[index]

# FCOALFF --------------------------------------------------------------------
BothMode_FCOALFF_SSP370 = xr.zeros_like(BC_COALFF_SSP370)
BC_COALFF_SSP370_copy = BC_COALFF_SSP370.copy()
for index in np.ndindex(BC_COALFF_SSP370.shape):
    if index in masked_indices_COALFF:
        BothMode_FCOALFF_SSP370[index] = BC_COALFF_SSP370_copy[index] * avg_iron_frac_valid
    else:
            BothMode_FCOALFF_SSP370[index] = (BothMode_FCOALFF_PD[index] / BC_COALFF_PD[index]) * BC_COALFF_SSP370_copy[index]

# FOIL --------------------------------------------------------------------
ocnfrac_values = ocnfrac.values  # Extract the raw values as a numpy array
ocnfrac_df = pd.DataFrame(ocnfrac_values, columns=ocnfrac['lon'].values, index=ocnfrac['lat'].values)
indices_gte_05 = []

for lat_idx in range(ocnfrac_df.shape[0]):  
    for lon_idx in range(ocnfrac_df.shape[1]):  
        if ocnfrac_df.iloc[lat_idx, lon_idx] >= 0.5:
            indices_gte_05.append((lat_idx, lon_idx))

count_gte_05 = len(indices_gte_05) # indices where ocnfrac was >= 0.5

BothMode_FOIL_SSP370 = xr.zeros_like(BC_OIL_SSP370)
BC_OIL_SSP370_copy = BC_OIL_SSP370.copy()
for index in np.ndindex(BC_OIL_SSP370.shape):
    if index in masked_indices_OIL:
   # For masked indices, check if the index is an ocean cell
        if index in indices_gte_05:
            # Multiply by avg_iron_frac_ocean if the index is in indices_gte_05
            BothMode_FOIL_SSP370[index] = BC_OIL_SSP370_copy[index] * avg_iron_frac_ocean
        else:
            # Otherwise, multiply by avg_iron_frac_valid
            BothMode_FOIL_SSP370[index] = BC_OIL_SSP370_copy[index] * avg_iron_frac_valid
    else:
        # For all other indices, use the Fe SPEWFUELS data to determine regional Fe emission frac 
            BothMode_FOIL_SSP370[index] = (BothMode_FOIL_PD[index] / BC_OIL_PD[index]) * BC_OIL_SSP370_copy[index]

# FWOOD --------------------------------------------------------------------
BothMode_FWOOD_SSP370 = xr.zeros_like(BC_WOOD_SSP370)
BC_WOOD_SSP370_copy = BC_WOOD_SSP370.copy()
for index in np.ndindex(BC_WOOD_SSP370.shape):
    if index in masked_indices_WOOD:
        BothMode_FWOOD_SSP370[index] = BC_WOOD_SSP370_copy[index] * avg_iron_frac_valid
    else:
            BothMode_FWOOD_SSP370[index] = (BothMode_FWOOD_PD[index] / BC_WOOD_PD[index]) * BC_WOOD_SSP370_copy[index]
# -------------------------------------------------------------------------------------------------------

# for manually checking dataframe and cell entries 
BothMode_FCOALBF_SSP370_df = pd.DataFrame(BothMode_FCOALBF_SSP370)
from itables import show
# Show the DataFrame interactively to check cells 
show(BothMode_FCOALBF_SSP370_df, scrollY=True, scrollX=True, maxRows=100, maxColumns=200) 

# ------- Size distribution (ratio of Coa to Fine mode Fe) is distinct for all sectors -----------------------
# -------- calculating the fine and coarse fractionation from the PD run to apply percentage each mode to the future data 
# -------- emissions are 0 over the ocean so the data over the ocean becomes NULL without changing to 0's 
CoaFrac_FCOALBF = CoaMed_FCOALBF_PD/(CoaMed_FCOALBF_PD + FineMed_FCOALBF_PD)
FineFrac_FCOALBF = FineMed_FCOALBF_PD/(CoaMed_FCOALBF_PD + FineMed_FCOALBF_PD)

CoaFrac_FCOALFF = CoaMed_FCOALFF_PD/(CoaMed_FCOALFF_PD + FineMed_FCOALFF_PD)
FineFrac_FCOALFF = FineMed_FCOALFF_PD/(CoaMed_FCOALFF_PD + FineMed_FCOALFF_PD)

CoaFrac_FOIL = CoaMed_FOIL_PD/(CoaMed_FOIL_PD + FineMed_FOIL_PD)
FineFrac_FOIL = FineMed_FOIL_PD/(CoaMed_FOIL_PD + FineMed_FOIL_PD)

CoaFrac_FWOOD = CoaMed_FWOOD_PD/(CoaMed_FWOOD_PD + FineMed_FWOOD_PD)
FineFrac_FWOOD = FineMed_FWOOD_PD/(CoaMed_FWOOD_PD + FineMed_FWOOD_PD)

average_coafrac_FCOALBF = CoaFrac_FCOALBF.mean()
print("avg coafrac_FCOALBF", average_coafrac_FCOALBF)
average_Finefrac_FCOALBF = FineFrac_FCOALBF.mean()
print("avg finefrac_FCOALBF", average_Finefrac_FCOALBF)
total_frac_FCOALBF = average_coafrac_FCOALBF + average_Finefrac_FCOALBF
print("avg_total_FCOALBF", total_frac_FCOALBF)

average_coafrac_FCOALFF = CoaFrac_FCOALFF.mean()
print("avg coafrac_FCOALFF", average_coafrac_FCOALFF)
average_Finefrac_FCOALFF = FineFrac_FCOALFF.mean()
print("avg finefrac_FCOALFF", average_Finefrac_FCOALFF)
total_frac_FCOALFF = average_coafrac_FCOALFF + average_Finefrac_FCOALFF
print("avg_total_FCOALFF", total_frac_FCOALFF)

average_coafrac_FOIL = CoaFrac_FOIL.mean()
print("avg coafrac_FOIL", average_coafrac_FOIL)
average_Finefrac_FOIL = FineFrac_FOIL.mean()
print("avg finefrac_FOIL", average_Finefrac_FOIL)
total_frac_FOIL = average_coafrac_FOIL + average_Finefrac_FOIL
print("avg_total_FOIL", total_frac_FOIL)

average_coafrac_FWOOD = CoaFrac_FWOOD.mean()
print("avg coafrac_FWOOD", average_coafrac_FWOOD)
average_Finefrac_FWOOD = FineFrac_FWOOD.mean()
print("avg finefrac_FWOOD", average_Finefrac_FWOOD)
total_frac_FWOOD = average_coafrac_FWOOD + average_Finefrac_FWOOD
print("avg_total_FWOOD", total_frac_FWOOD)
# totals all add to one, passed this test

# now to size fractionate the FU Fe, using individual grid cells where possible, but average size fractionation by sector where data is missing
# BIOFUELS ------------------------------------------------------------------------------------------------
CoaMed_FCOALBF_SSP370 = xr.zeros_like(BothMode_FCOALBF_SSP370)
for index in np.ndindex(BothMode_FCOALBF_SSP370.shape):
    if index in masked_indices_COALBF:
        # For masked indices where data does not exist for PD files, multiply the cell by the avg aerosol mode fraction for that sector
        CoaMed_FCOALBF_SSP370[index] = BothMode_FCOALBF_SSP370[index] * average_coafrac_FCOALBF
    else:
        # For all other indices, use the Fe SPEWFUELS data to determine regional Fe emission frac 
        CoaMed_FCOALBF_SSP370[index] = BothMode_FCOALBF_SSP370[index] * CoaFrac_FCOALBF[index]

FineMed_FCOALBF_SSP370 = xr.zeros_like(BothMode_FCOALBF_SSP370)
for index in np.ndindex(BothMode_FCOALBF_SSP370.shape):
    if index in masked_indices_COALBF:
        # For masked indices where data does not exist for PD files, multiply the cell by the avg aerosol mode fraction for that sector
        FineMed_FCOALBF_SSP370[index] = BothMode_FCOALBF_SSP370[index] * average_Finefrac_FCOALBF
    else:
        # For all other indices, use the Fe SPEWFUELS data to determine regional Fe emission frac 
        FineMed_FCOALBF_SSP370[index] = BothMode_FCOALBF_SSP370[index] * FineFrac_FCOALBF[index]

# FOSSIL FUELS ------------------------------------------------------------------------------------------------
CoaMed_FCOALFF_SSP370 = xr.zeros_like(BothMode_FCOALFF_SSP370)
for index in np.ndindex(BothMode_FCOALFF_SSP370.shape):
    if index in masked_indices_COALBF:
        # For masked indices where data does not exist for PD files, multiply the cell by the avg aerosol mode fraction for that sector
        CoaMed_FCOALFF_SSP370[index] = BothMode_FCOALFF_SSP370[index] * average_coafrac_FCOALFF
    else:
        # For all other indices, use the Fe SPEWFUELS data to determine regional Fe emission frac 
        CoaMed_FCOALFF_SSP370[index] = BothMode_FCOALFF_SSP370[index] * CoaFrac_FCOALFF[index]

FineMed_FCOALFF_SSP370 = xr.zeros_like(BothMode_FCOALFF_SSP370)
for index in np.ndindex(BothMode_FCOALFF_SSP370.shape):
    if index in masked_indices_COALBF:
        # For masked indices where data does not exist for PD files, multiply the cell by the avg aerosol mode fraction for that sector
        FineMed_FCOALFF_SSP370[index] = BothMode_FCOALFF_SSP370[index] * average_Finefrac_FCOALFF
    else:
        # For all other indices, use the Fe SPEWFUELS data to determine regional Fe emission frac 
        FineMed_FCOALFF_SSP370[index] = BothMode_FCOALFF_SSP370[index] * FineFrac_FCOALFF[index]

# OIL ------------------------------------------------------------------------------------------------
CoaMed_FOIL_SSP370 = xr.zeros_like(BothMode_FOIL_SSP370)
for index in np.ndindex(BothMode_FOIL_SSP370.shape):
    if index in masked_indices_COALBF:
        # For masked indices where data does not exist for PD files, multiply the cell by the avg aerosol mode fraction for that sector
        CoaMed_FOIL_SSP370[index] = BothMode_FOIL_SSP370[index] * average_coafrac_FOIL
    else:
        # For all other indices, use the Fe SPEWFUELS data to determine regional Fe emission frac 
        CoaMed_FOIL_SSP370[index] = BothMode_FOIL_SSP370[index] * CoaFrac_FOIL[index]

FineMed_FOIL_SSP370 = xr.zeros_like(BothMode_FOIL_SSP370)
for index in np.ndindex(BothMode_FOIL_SSP370.shape):
    if index in masked_indices_COALBF:
        # For masked indices where data does not exist for PD files, multiply the cell by the avg aerosol mode fraction for that sector
        FineMed_FOIL_SSP370[index] = BothMode_FOIL_SSP370[index] * average_Finefrac_FOIL
    else:
        # For all other indices, use the Fe SPEWFUELS data to determine regional Fe emission frac 
        FineMed_FOIL_SSP370[index] = BothMode_FOIL_SSP370[index] * FineFrac_FOIL[index]

# WOOD ------------------------------------------------------------------------------------------------
CoaMed_FWOOD_SSP370 = xr.zeros_like(BothMode_FWOOD_SSP370)
for index in np.ndindex(BothMode_FWOOD_SSP370.shape):
    if index in masked_indices_COALBF:
        # For masked indices where data does not exist for PD files, multiply the cell by the avg aerosol mode fraction for that sector
        CoaMed_FWOOD_SSP370[index] = BothMode_FWOOD_SSP370[index] * average_coafrac_FWOOD
    else:
        # For all other indices, use the Fe SPEWFUELS data to determine regional Fe emission frac 
        CoaMed_FWOOD_SSP370[index] = BothMode_FWOOD_SSP370[index] * CoaFrac_FWOOD[index]

FineMed_FWOOD_SSP370 = xr.zeros_like(BothMode_FWOOD_SSP370)
for index in np.ndindex(BothMode_FWOOD_SSP370.shape):
    if index in masked_indices_COALBF:
        # For masked indices where data does not exist for PD files, multiply the cell by the avg aerosol mode fraction for that sector
        FineMed_FWOOD_SSP370[index] = BothMode_FWOOD_SSP370[index] * average_Finefrac_FWOOD
    else:
        # For all other indices, use the Fe SPEWFUELS data to determine regional Fe emission frac 
        FineMed_FWOOD_SSP370[index] = BothMode_FWOOD_SSP370[index] * FineFrac_FWOOD[index]

# for manually checking dataframe and cell entries 
#FineMed_FCOALBF_SSP370_df = pd.DataFrame(FineMed_FCOALBF_SSP370)
#from itables import show
# Show the DataFrame interactively to check cells 
#show(FineMed_FCOALBF_SSP370_df, scrollY=True, scrollX=True, maxRows=100, maxColumns=200) 

# REPEATING FRACTIONATION FOR BLACK CARBON SSP370
# BIOFUELS ------------------------------------------------------------------------------------------------
CoaMed_BCCOALBF_SSP370 = xr.zeros_like(BC_COALBF_PD)
for index in np.ndindex(BC_COALBF_PD.shape):
    if index in masked_indices_COALBF:
# # For masked indices where data does not exist for PD files, multiply the cell by the avg aerosol mode fraction for that sector
        CoaMed_BCCOALBF_SSP370[index] = BC_COALBF_PD[index] * average_coafrac_FCOALBF
    else:
#        # For all other indices, use the Fe SPEWFUELS data to determine regional Fe emission frac 
        CoaMed_BCCOALBF_SSP370[index] = BC_COALBF_PD[index] * CoaFrac_FCOALBF[index]

FineMed_BCCOALBF_SSP370 = xr.zeros_like(BothMode_FCOALBF_SSP370)
for index in np.ndindex(BothMode_FCOALBF_SSP370.shape):
    if index in masked_indices_COALBF:
 #       # For masked indices where data does not exist for PD files, multiply the cell by the avg aerosol mode fraction for that sector
        FineMed_BCCOALBF_SSP370[index] = BC_COALBF_PD[index] * average_Finefrac_FCOALBF
    else:
 #       # For all other indices, use the Fe SPEWFUELS data to determine regional Fe emission frac 
        FineMed_BCCOALBF_SSP370[index] = BC_COALBF_PD[index] * FineFrac_FCOALBF[index]

# FOSSIL FUELS ------------------------------------------------------------------------------------------------
CoaMed_BCCOALFF_SSP370 = xr.zeros_like(BC_COALFF_PD)
for index in np.ndindex(BC_COALFF_PD.shape):
    if index in masked_indices_COALBF:
#        # For masked indices where data does not exist for PD files, multiply the cell by the avg aerosol mode fraction for that sector
        CoaMed_BCCOALFF_SSP370[index] = BC_COALFF_PD[index] * average_coafrac_FCOALFF
    else:
 #       # For all other indices, use the Fe SPEWFUELS data to determine regional Fe emission frac 
        CoaMed_BCCOALFF_SSP370[index] = BC_COALFF_PD[index] * CoaFrac_FCOALFF[index]

FineMed_BCCOALFF_SSP370 = xr.zeros_like(BothMode_FCOALFF_SSP370)
for index in np.ndindex(BothMode_FCOALFF_SSP370.shape):
    if index in masked_indices_COALBF:
 #       # For masked indices where data does not exist for PD files, multiply the cell by the avg aerosol mode fraction for that sector
        FineMed_BCCOALFF_SSP370[index] = BC_COALFF_PD[index] * average_Finefrac_FCOALFF
    else:
        # For all other indices, use the Fe SPEWFUELS data to determine regional Fe emission frac 
        FineMed_BCCOALFF_SSP370[index] = BC_COALFF_PD[index] * FineFrac_FCOALFF[index]

# OIL ------------------------------------------------------------------------------------------------
CoaMed_BCOIL_SSP370 = xr.zeros_like(BC_OIL_PD)
for index in np.ndindex(BC_OIL_PD.shape):
    if index in masked_indices_COALBF:
        # For masked indices where data does not exist for PD files, multiply the cell by the avg aerosol mode fraction for that sector
        CoaMed_BCOIL_SSP370[index] = BC_OIL_PD[index] * average_coafrac_FOIL
    else:
   #     # For all other indices, use the Fe SPEWFUELS data to determine regional Fe emission frac 
        CoaMed_BCOIL_SSP370[index] = BC_OIL_PD[index] * CoaFrac_FOIL[index]

FineMed_BCOIL_SSP370 = xr.zeros_like(BC_OIL_PD)
for index in np.ndindex(BC_OIL_PD.shape):
    if index in masked_indices_COALBF:
 #       # For masked indices where data does not exist for PD files, multiply the cell by the avg aerosol mode fraction for that sector
        FineMed_BCOIL_SSP370[index] = BC_OIL_PD[index] * average_Finefrac_FOIL
    else:
        # For all other indices, use the Fe SPEWFUELS data to determine regional Fe emission frac 
        FineMed_BCOIL_SSP370[index] = BC_OIL_PD[index] * FineFrac_FOIL[index]

# WOOD ------------------------------------------------------------------------------------------------
CoaMed_BCWOOD_SSP370 = xr.zeros_like(BC_WOOD_PD)
for index in np.ndindex(BC_WOOD_PD.shape):
    if index in masked_indices_COALBF:
        # For masked indices where data does not exist for PD files, multiply the cell by the avg aerosol mode fraction for that sector
        CoaMed_BCWOOD_SSP370[index] = BC_WOOD_PD[index] * average_coafrac_FWOOD
    else:
        # For all other indices, use the Fe SPEWFUELS data to determine regional Fe emission frac 
        CoaMed_BCWOOD_SSP370[index] = BC_WOOD_PD[index] * CoaFrac_FWOOD[index]

FineMed_BCWOOD_SSP370 = xr.zeros_like(BC_WOOD_PD)
for index in np.ndindex(BC_WOOD_PD.shape):
    if index in masked_indices_COALBF:
        # For masked indices where data does not exist for PD files, multiply the cell by the avg aerosol mode fraction for that sector
        FineMed_BCWOOD_SSP370[index] = BC_WOOD_PD[index] * average_Finefrac_FWOOD
    else:
        # For all other indices, use the Fe SPEWFUELS data to determine regional Fe emission frac 
        FineMed_BCWOOD_SSP370[index] = BC_WOOD_PD[index] * FineFrac_FWOOD[index]

# ----- compared these values to Rathod et al 2019 Figure 2 and the ratios Coa:Fine are approximately correct for each sector, they passed this test -------
# Plots testing fractions 
FineMed_FWOOD_SSP370.plot(subplot_kws={'projection': ccrs.PlateCarree()},transform=ccrs.PlateCarree(), cbar_kwargs={'label': 'Fraction'})
plt.title("FineFrac_FCOALBF")
plt.gca().coastlines()
plt.show() 

FineMed_BCWOOD_SSP370.plot(subplot_kws={'projection': ccrs.PlateCarree()},transform=ccrs.PlateCarree(), cbar_kwargs={'label': 'Fraction'})
#plt.title("CoaFrac_FCOALBF")
#plt.gca().coastlines()
#plt.show()

# --------- separating BC emissions into fine and coarse fraction based on specific grid cell ratios from Fe file 
# ------ need to ensure dimensions are the same or this wont work (no straggling sectors or time coords) ------------
# ----- BIOFUELS

CoaMed_BCCOALBF_SSP370 = CoaFrac_FCOALBF * BC_COALBF_SSP370
#CoaMed_BCCOALBF_SSP370 = CoaMed_BCCOALBF_SSP370.fillna(0.0)
#CoaMed_BCCOALBF_SSP370 = CoaMed_BCCOALBF_SSP370.where(np.isfinite(CoaMed_BCCOALBF_SSP370), 0.0)

FineMed_BCCOALBF_SSP370 = FineFrac_FCOALBF * BC_COALBF_SSP370 
#FineMed_BCCOALBF_SSP370 = FineMed_BCCOALBF_SSP370.fillna(0.0)
#FineMed_BCCOALBF_SSP370 = FineMed_BCCOALBF_SSP370.where(np.isfinite(FineMed_BCCOALBF_SSP370), 0.0)

#total_BCCOALBF_SSP370 = CoaMed_BCCOALBF_SSP370 + FineMed_BCCOALBF_SSP370

CoaMed_BCCOALBF_PD = CoaFrac_FCOALBF * BC_COALBF_PD
CoaMed_BCCOALBF_PD = CoaMed_BCCOALBF_PD.fillna(0.0)
CoaMed_BCCOALBF_PD = CoaMed_BCCOALBF_PD.where(np.isfinite(CoaMed_BCCOALBF_PD), 0.0)

FineMed_BCCOALBF_PD = FineFrac_FCOALBF * BC_COALBF_PD
FineMed_BCCOALBF_PD = FineMed_BCCOALBF_PD.fillna(0.0)
FineMed_BCCOALBF_PD = FineMed_BCCOALBF_PD.where(np.isfinite(FineMed_BCCOALBF_PD), 0.0)

#total_BCCOALBF_PD = CoaMed_BCCOALBF_PD + FineMed_BCCOALBF_PD

#print("OG_BC", BC_COALBF_PD)
#print("CoaFrac", CoaFrac_FCOALBF)
#print("CoaMed_BCCOALBF_SSP370", CoaMed_BCCOALBF_SSP370) 
# ----- they do add to total QA/QC -------

#print("CoaMed_BC", CoaMed_BCCOALBF_PD.sum())
#print("FineMed_BC", FineMed_BCCOALBF_PD.sum())
#print("Total_BC", total_BCCOALBF_PD.sum()) 
# ----- they do add to total QA/QC -------

# ----- wrote out to test if BC >> Fe, confirmed for at least BFCOAL -----------
# CoaMed_BCCOALBF_PD.to_netcdf("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\CoaMed_BCCOALBF_PD.nc")

# ----- FOSSIL FUELS
CoaMed_BCCOALFF_SSP370 = CoaFrac_FCOALFF * BC_COALFF_SSP370 
#CoaMed_BCCOALFF_SSP370 = CoaMed_BCCOALFF_SSP370.fillna(0.0)
#CoaMed_BCCOALFF_SSP370 = CoaMed_BCCOALFF_SSP370.where(np.isfinite(CoaMed_BCCOALFF_SSP370), 0.0)

FineMed_BCCOALFF_SSP370 = FineFrac_FCOALFF * BC_COALFF_SSP370 
#FineMed_BCCOALFF_SSP370 = FineMed_BCCOALFF_SSP370.fillna(0.0)
#FineMed_BCCOALFF_SSP370 = FineMed_BCCOALFF_SSP370.where(np.isfinite(FineMed_BCCOALFF_SSP370), 0.0)

#total_BCCOALFF_SSP370 = CoaMed_BCCOALFF_SSP370  + FineMed_BCCOALFF_SSP370 

CoaMed_BCCOALFF_PD = CoaFrac_FCOALFF * BC_COALFF_PD 
CoaMed_BCCOALFF_PD = CoaMed_BCCOALFF_PD.fillna(0.0)
CoaMed_BCCOALFF_PD  = CoaMed_BCCOALFF_PD.where(np.isfinite(CoaMed_BCCOALFF_PD), 0.0)

FineMed_BCCOALFF_PD = FineFrac_FCOALFF * BC_COALFF_PD 
FineMed_BCCOALFF_PD = FineMed_BCCOALFF_PD.fillna(0.0)
FineMed_BCCOALFF_PD = FineMed_BCCOALFF_PD.where(np.isfinite(FineMed_BCCOALFF_PD), 0.0)

#total_BCCOALFF_PD = CoaMed_BCCOALFF_PD + FineMed_BCCOALFF_PD

# ----- OIL
CoaMed_BCOIL_SSP370 = CoaFrac_FOIL * BC_OIL_SSP370
#CoaMed_BCOIL_SSP370 = CoaMed_BCOIL_SSP370.fillna(0.0)
#CoaMed_BCOIL_SSP370 = CoaMed_BCOIL_SSP370.where(np.isfinite(CoaMed_BCOIL_SSP370), 0.0)

FineMed_BCOIL_SSP370 = FineFrac_FOIL* BC_OIL_SSP370
#FineMed_BCOIL_SSP370 = FineMed_BCOIL_SSP370.fillna(0.0)
#FineMed_BCOIL_SSP370 = FineMed_BCOIL_SSP370.where(np.isfinite(FineMed_BCOIL_SSP370), 0.0)

#total_BCOIL_SSP370 = CoaMed_BCOIL_SSP370 + FineMed_BCOIL_SSP370

CoaMed_BCOIL_PD = CoaFrac_FOIL * BC_OIL_PD
CoaMed_BCOIL_PD = CoaMed_BCOIL_PD.fillna(0.0)
CoaMed_BCOIL_PD = CoaMed_BCOIL_PD.where(np.isfinite(CoaMed_BCOIL_PD), 0.0)

FineMed_BCOIL_PD = FineFrac_FOIL * BC_OIL_PD
FineMed_BCOIL_PD = FineMed_BCOIL_PD.fillna(0.0)
FineMed_BCOIL_PD = FineMed_BCOIL_PD.where(np.isfinite(FineMed_BCOIL_PD), 0.0)

#total_BCOIL_PD = CoaMed_BCOIL_PD + FineMed_BCOIL_PD

# ----- WOOD
CoaMed_BCWOOD_SSP370 = CoaFrac_FWOOD * BC_WOOD_SSP370
#CoaMed_BCWOOD_SSP370 = CoaMed_BCWOOD_SSP370.fillna(0.0)
#CoaMed_BCWOOD_SSP370 = CoaMed_BCWOOD_SSP370.where(np.isfinite(CoaMed_BCWOOD_SSP370), 0.0)

FineMed_BCWOOD_SSP370 = FineFrac_FWOOD * BC_WOOD_SSP370
#FineMed_BCWOOD_SSP370 = FineMed_BCWOOD_SSP370.fillna(0.0)
#FineMed_BCWOOD_SSP370 = FineMed_BCWOOD_SSP370.where(np.isfinite(CoaMed_BCWOOD_SSP370), 0.0)

#total_BCWOOD_SSP370 = CoaMed_BCWOOD_SSP370 + FineMed_BCWOOD_SSP370

CoaMed_BCWOOD_PD = CoaFrac_FWOOD * BC_WOOD_PD
CoaMed_BCWOOD_PD = CoaMed_BCWOOD_PD.fillna(0.0)
CoaMed_BCWOOD_PD = CoaMed_BCWOOD_PD.where(np.isfinite(CoaMed_BCWOOD_PD), 0.0)

FineMed_BCWOOD_PD = FineFrac_FWOOD * BC_WOOD_PD
FineMed_BCWOOD_PD = FineMed_BCWOOD_PD.fillna(0.0)
FineMed_BCWOOD_PD = FineMed_BCWOOD_PD.where(np.isfinite(FineMed_BCWOOD_PD), 0.0)

#total_BCWOOD_PD = CoaMed_BCWOOD_PD + FineMed_BCWOOD_PD

# ran this for every combination to check if fine frac + coa frac added to total QA/QC passed for each sector
# if mean error was <e-14 assume rounding differences
# print("OG", BC_WOOD_SSP370.mean().item())
# print("Coa", CoaMed_BCWOOD_SSP370.mean().item())
# print("Fine", FineMed_BCWOOD_SSP370.mean().item())
# print("added", total_BCWOOD_SSP370.mean().item()) 
# difference = BC_WOOD_SSP370 - total_BCWOOD_SSP370
# print("Mean difference:", difference.mean().item()) 
# results for differences below
# COALBF_SSP370 = e-15
# COALBF_PD = e-15
# COALFF_SSP370 = e-15
# COALFF_PD = e-15
# OIL_SSP370 = e-16
# OIL_PD = e-16
# WOOD_SSP370 = e-15
# WOOD_PD = e-15

# Checking ratios of Fe to BC FOR PRESENT DAY 
# Calculating mean Fe to BC ratio for PD data -- to ensure Fe << BC
#FeBC_ratio = CoaMed_FWOOD_PD/CoaMed_BCWOOD_PD
#FeBC_ratio = FeBC_ratio.fillna(0.0)
#FeBC_ratio = FeBC_ratio.where(np.isfinite(FeBC_ratio), 0.0)
#nan_count_1 = FeBC_ratio.isnull().sum().values
#nan_count_2 = CoaMed_BCWOOD_PD.isnull().sum().values
#print(nan_count_1)
#print(nan_count_2)
#valid_values = FeBC_ratio.where(np.isfinite(FeBC_ratio), drop=True) # all this does is replace infinities with NaNs

#FeBC_ratio_nanmask = FeBC_ratio.isnull()

# Plotting regions where Fe > BC -- will present to group meeting first week of January 2025
#exceedance_indices = []
#xarr1 = CoaMed_FWOOD_PD.where(np.isfinite(CoaMed_FWOOD_PD))
#xarr2 = CoaMed_BCWOOD_PD.where(np.isfinite(CoaMed_BCWOOD_PD))

# used this to do a quick check to make sure different values were coming up uniquely for each sector
#from itables import show
#xarr1_df  = pd.DataFrame(CoaMed_BCWOOD_PD)
# Show the DataFrame interactively -- to see what inf and NaNs are doing
# show(xarr1_df, scrollY=True, scrollX=True, maxRows=100, maxColumns=200) # looks good, data is dif between sectors

#for i in range(xarr1.shape[0]):
#   for j in range(xarr1.shape[1]):
#        if xarr1[i, j] >= xarr2[i, j]:
#           exceedance_indices.append((i, j))

#valid_mask = np.isfinite(xarr1) & np.isfinite(xarr2)
#exceedance_map = ((xarr1 > xarr2) & valid_mask).astype(int)
#exceedance_percentage = (exceedance_map.sum() / valid_mask.sum()) * 100
#print(f"Percentage of Exceedance: {exceedance_percentage:.2f}%")
#NaN_percentage = (nan_count_xarr2/55296) * 100
#print(f"Percentage of NaNs: {NaN_percentage:.2f}%")

#exceedance_map = (FeBC_ratio > 1).astype(int)

#lat = FeBC_ratio["lat"]
#lon = FeBC_ratio["lon"]

#import cartopy.feature as cfeature
#import matplotlib.colors as mcolors
#custom_cmap = mcolors.ListedColormap(["white", "blue"])
#fig, ax = plt.subplots(
#    figsize=(10, 8),
#    subplot_kw={"projection": ccrs.Robinson()},)

#im = ax.pcolormesh(
#   lon, lat, exceedance_map.values,
#   cmap=custom_cmap,  # Flipped colormap
 #   shading="auto",
 #  transform=ccrs.PlateCarree())

#ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
#ax.add_feature(cfeature.BORDERS, linestyle=":")
#gl = ax.gridlines(draw_labels=True, linestyle="--", color="gray", alpha=0.5)
#gl.top_labels = False
#gl.right_labels = False
#ax.set_title("grid cells where COARSE WOOD Fe exceeds BC (Present Day)", fontsize=14)
#plt.show()

# Count the number of NaN values
#nan_count = FeBC_ratio.isnull().sum().values
# Mask the valid values that are between 0 and 1
#ratios_in_range = (valid_values <= 1)
# Count the number of valid values in the range 0-1
#count_in_range = ratios_in_range.sum().values  # Extract scalar value
# Count the number of valid values outside the range [0, 1]
#other_count = valid_values.count().values - count_in_range

#print(f"% NaNs: {nan_count/FeBC_ratio.size}")
#print(f"% Fe does not exceed BC: {count_in_range/FeBC_ratio.size}")
#print(f"% Fe exceeds BC: {other_count/FeBC_ratio.size}")

# after extensive testing, it appears that any discrepancies between calculating Fe/BC and Fe>BC is due to how NaNs are being handled, so I should go ahead and move forward to calculate budgets.

# data completeness seems low -- all the areas where BC is zero when looking at plots -- replace NANs with 0, 
# but it appears to match the original BC files per emissions mostly from terrestrial sites which are only 30% of the Earth 
# check for other sectors -- looks good 

# This is where I need to double check my work, looking at 
# let's break this down into two steps -- find ratio of Fe:BC PD then multiply through by BC SSP370
# Calculating Fe_SSP370 emission values from ratio of FePD/BCPD to BCFU
# FCOALBF --------------------------------------------------------------------
#CoaMed_FCOALBF_SSP370 = (CoaMed_FCOALBF_PD/CoaMed_BCCOALBF_PD) * CoaMed_BCCOALBF_SSP370
#CoaMed_FCOALBF_SSP370 = CoaMed_FCOALBF_SSP370.fillna(0.0) # don't replace NANs until very end 
#FineMed_FCOALBF_SSP370 = (FineMed_FCOALBF_PD/FineMed_BCCOALBF_PD) * FineMed_BCCOALBF_SSP370
#FineMed_FCOALBF_SSP370 = FineMed_FCOALBF_SSP370.fillna(0.0)

# Find ratio (fraction) of Fe:BC PD
#CoaMed_COALBF_FeBC_frac_PD = (CoaMed_FCOALBF_PD/CoaMed_BCCOALBF_PD)
# replacing NaNs and Infinities is absolutely imperative for the ratio of Fe:BC values (at least averages checked in panoply) to be conserved between PD and SSP370
#CoaMed_COALBF_FeBC_frac_PD = CoaMed_COALBF_FeBC_frac_PD.fillna(0.0)
#CoaMed_COALBF_FeBC_frac_PD = CoaMed_COALBF_FeBC_frac_PD.where(np.isfinite(CoaMed_COALBF_FeBC_frac_PD), 0.0)

#FineMed_COALBF_FeBC_frac_PD = (FineMed_FCOALBF_PD/FineMed_BCCOALBF_PD)
#FineMed_COALBF_FeBC_frac_PD = FineMed_COALBF_FeBC_frac_PD.fillna(0.0)
#FineMed_COALBF_FeBC_frac_PD = FineMed_COALBF_FeBC_frac_PD.where(np.isfinite(FineMed_COALBF_FeBC_frac_PD), 0.0)

#CoaMed_FCOALBF_SSP370 = CoaMed_COALBF_FeBC_frac_PD * CoaMed_BCCOALBF_SSP370
#CoaMed_FCOALBF_SSP370 = CoaMed_FCOALBF_SSP370.fillna(0.0)
#CoaMed_FCOALBF_SSP370 = CoaMed_FCOALBF_SSP370.where(np.isfinite(CoaMed_FCOALBF_SSP370), 0.0)

#FineMed_FCOALBF_SSP370 = FineMed_COALBF_FeBC_frac_PD * FineMed_BCCOALBF_SSP370
#FineMed_FCOALBF_SSP370 = FineMed_FCOALBF_SSP370.fillna(0.0)
#FineMed_FCOALBF_SSP370 = FineMed_FCOALBF_SSP370.where(np.isfinite(FineMed_FCOALBF_SSP370), 0.0)

#CoaMed_COALBF_FeBC_frac_SSP370 = (CoaMed_FCOALBF_SSP370/CoaMed_BCCOALBF_SSP370)
#CoaMed_COALBF_FeBC_frac_SSP370 = CoaMed_COALBF_FeBC_frac_SSP370.fillna(0.0)
#CoaMed_COALBF_FeBC_frac_SSP370 = CoaMed_COALBF_FeBC_frac_SSP370.where(np.isfinite(CoaMed_COALBF_FeBC_frac_SSP370), 0.0)

#FineMed_COALBF_FeBC_frac_SSP370 = (FineMed_FCOALBF_SSP370/FineMed_BCCOALBF_SSP370)
#FineMed_COALBF_FeBC_frac_SSP370 = FineMed_COALBF_FeBC_frac_SSP370.fillna(0.0)
#FineMed_COALBF_FeBC_frac_SSP370 = FineMed_COALBF_FeBC_frac_SSP370.where(np.isfinite(FineMed_COALBF_FeBC_frac_SSP370), 0.0)

#Check to see if ratios look the same -- there were some subtle differences when infinities weren't properly removed as discovered in Panoply

#CoaMed_COALBF_FeBC_frac_PD.to_netcdf("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\FeBC_ratio_PD.nc")
#CoaMed_COALBF_FeBC_frac_SSP370.to_netcdf("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\FeBC_ratio_SSP370.nc")
# there appear to be many infinities and NaNs (that vary between PD and SSP370 -- I need to ensure that NaNs and Infinities are BOTH recoded as zero for emissions)

# FCOALFF --------------------------------------------------------------------
#CoaMed_FCOALFF_SSP370 = (CoaMed_FCOALFF_PD/CoaMed_BCCOALFF_PD) * CoaMed_BCCOALFF_SSP370
#CoaMed_FCOALFF_SSP370 = CoaMed_FCOALFF_SSP370.fillna(0.0)
#FineMed_FCOALFF_SSP370 = (FineMed_FCOALFF_PD/FineMed_BCCOALFF_PD) * FineMed_BCCOALFF_SSP370
#FineMed_FCOALFF_SSP370 = FineMed_FCOALFF_SSP370.fillna(0.0)

#CoaMed_COALFF_FeBC_frac_PD = (CoaMed_FCOALFF_PD/CoaMed_BCCOALFF_PD)
#CoaMed_COALFF_FeBC_frac_PD = CoaMed_COALFF_FeBC_frac_PD.fillna(0.0)
#CoaMed_COALFF_FeBC_frac_PD = CoaMed_COALFF_FeBC_frac_PD.where(np.isfinite(CoaMed_COALFF_FeBC_frac_PD), 0.0)

#FineMed_COALFF_FeBC_frac_PD = (FineMed_FCOALFF_PD/FineMed_BCCOALFF_PD)
#FineMed_COALFF_FeBC_frac_PD = FineMed_COALFF_FeBC_frac_PD.fillna(0.0)
#FineMed_COALFF_FeBC_frac_PD = FineMed_COALFF_FeBC_frac_PD.where(np.isfinite(FineMed_COALFF_FeBC_frac_PD), 0.0)

#CoaMed_FCOALFF_SSP370 = CoaMed_COALFF_FeBC_frac_PD * CoaMed_BCCOALFF_SSP370
#CoaMed_FCOALFF_SSP370 = CoaMed_FCOALFF_SSP370.fillna(0.0)
#CoaMed_FCOALFF_SSP370 = CoaMed_FCOALFF_SSP370.where(np.isfinite(FineMed_FCOALBF_SSP370), 0.0)

#FineMed_FCOALFF_SSP370 = FineMed_COALFF_FeBC_frac_PD * FineMed_BCCOALFF_SSP370
#FineMed_FCOALFF_SSP370 = FineMed_FCOALFF_SSP370.fillna(0.0)
#FineMed_FCOALFF_SSP370 = FineMed_FCOALFF_SSP370.where(np.isfinite(FineMed_FCOALFF_SSP370), 0.0)

#CoaMed_COALFF_FeBC_frac_SSP370 = (CoaMed_FCOALFF_SSP370/CoaMed_BCCOALFF_SSP370)
#CoaMed_COALFF_FeBC_frac_SSP370 = CoaMed_COALFF_FeBC_frac_SSP370.fillna(0.0)
#CoaMed_COALFF_FeBC_frac_SSP370 = CoaMed_COALFF_FeBC_frac_SSP370.where(np.isfinite(CoaMed_COALFF_FeBC_frac_SSP370), 0.0)

#FineMed_COALFF_FeBC_frac_SSP370 = (FineMed_FCOALFF_SSP370/FineMed_BCCOALFF_SSP370)
#FineMed_COALFF_FeBC_frac_SSP370 = FineMed_COALFF_FeBC_frac_SSP370.fillna(0.0)
#FineMed_COALFF_FeBC_frac_SSP370 = FineMed_COALFF_FeBC_frac_SSP370.where(np.isfinite(FineMed_COALFF_FeBC_frac_SSP370), 0.0)

# these below are to test ratios of Fe to BC on map and examine in Panoply
#FineMed_COALFF_FeBC_frac_SSP370.to_netcdf("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\FeBC_ratio_PD.nc")
#FineMed_COALFF_FeBC_frac_PD.to_netcdf("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal #Fly Ash\\data\\FeBC_ratio_SSP370.nc")

# FOIL --------------------------------------------------------------------
#CoaMed_FOIL_SSP370 = (CoaMed_FOIL_PD/CoaMed_BCOIL_PD) * CoaMed_BCOIL_SSP370
#CoaMed_FOIL_SSP370 = CoaMed_FOIL_SSP370.fillna(0.0)
#FineMed_FOIL_SSP370 = (FineMed_FOIL_PD/FineMed_BCOIL_PD) * FineMed_BCOIL_SSP370
#FineMed_FOIL_SSP370 = FineMed_FOIL_SSP370.fillna(0.0)

#CoaMed_OIL_FeBC_frac_PD = (CoaMed_FOIL_PD/CoaMed_BCOIL_PD)
#CoaMed_OIL_FeBC_frac_PD = CoaMed_OIL_FeBC_frac_PD.fillna(0.0)
#CoaMed_OIL_FeBC_frac_PD = CoaMed_OIL_FeBC_frac_PD.where(np.isfinite(CoaMed_OIL_FeBC_frac_PD), 0.0)

#FineMed_OIL_FeBC_frac_PD = (FineMed_FOIL_PD/FineMed_BCOIL_PD)
#FineMed_OIL_FeBC_frac_PD = FineMed_OIL_FeBC_frac_PD.fillna(0.0)
#FineMed_OIL_FeBC_frac_PD = FineMed_OIL_FeBC_frac_PD.where(np.isfinite(FineMed_OIL_FeBC_frac_PD), 0.0)

#CoaMed_FOIL_SSP370 = CoaMed_OIL_FeBC_frac_PD * CoaMed_BCOIL_SSP370
#CoaMed_FOIL_SSP370 = CoaMed_FOIL_SSP370.fillna(0.0)
#CoaMed_FOIL_SSP370 = CoaMed_FOIL_SSP370.where(np.isfinite(CoaMed_FOIL_SSP370), 0.0)

#FineMed_FOIL_SSP370 = FineMed_OIL_FeBC_frac_PD * FineMed_BCOIL_SSP370
#FineMed_FOIL_SSP370 = FineMed_FOIL_SSP370.fillna(0.0)
#FineMed_FOIL_SSP370= FineMed_FOIL_SSP370.where(np.isfinite(FineMed_FOIL_SSP370), 0.0)

#CoaMed_OIL_FeBC_frac_SSP370 = (CoaMed_FOIL_SSP370/CoaMed_BCOIL_SSP370)
#CoaMed_OIL_FeBC_frac_SSP370 = CoaMed_OIL_FeBC_frac_SSP370.fillna(0.0)
#CoaMed_OIL_FeBC_frac_SSP370 = CoaMed_OIL_FeBC_frac_SSP370.where(np.isfinite(CoaMed_OIL_FeBC_frac_SSP370), 0.0)

#FineMed_OIL_FeBC_frac_SSP370 = (FineMed_FOIL_SSP370/FineMed_BCOIL_SSP370)
#FineMed_OIL_FeBC_frac_SSP370 = FineMed_OIL_FeBC_frac_SSP370.fillna(0.0)
#FineMed_OIL_FeBC_frac_SSP370 = FineMed_OIL_FeBC_frac_SSP370.where(np.isfinite(FineMed_OIL_FeBC_frac_SSP370), 0.0)

# WOOD -----------------------------------------------------------------------
#CoaMed_FWOOD_SSP370 = (CoaMed_FWOOD_PD/CoaMed_BCWOOD_PD) * CoaMed_BCWOOD_SSP370
#CoaMed_FWOOD_SSP370 = CoaMed_FWOOD_SSP370.fillna(0.0)
#FineMed_FWOOD_SSP370 = (FineMed_FWOOD_PD/FineMed_BCWOOD_PD) * FineMed_BCWOOD_SSP370
#FineMed_FWOOD_SSP370 = FineMed_FWOOD_SSP370.fillna(0.0)

#CoaMed_WOOD_FeBC_frac_PD = (CoaMed_FWOOD_PD/CoaMed_BCWOOD_PD)
# Do I get differences using this code from StackOverflow vs. the NaN and is infinite removals -- NO
#CoaMed_WOOD_FeBC_frac_PD = np.nan_to_num(CoaMed_WOOD_FeBC_frac_PD, nan=0.0, posinf=0.0, neginf=0.0)
#CoaMed_WOOD_FeBC_frac_PD = CoaMed_WOOD_FeBC_frac_PD.fillna(0.0)
#CoaMed_WOOD_FeBC_frac_PD = CoaMed_WOOD_FeBC_frac_PD.where(np.isfinite(CoaMed_WOOD_FeBC_frac_PD), 0.0)

#FineMed_WOOD_FeBC_frac_PD = (FineMed_FWOOD_PD/FineMed_BCWOOD_PD)
#FineMed_WOOD_FeBC_frac_PD = FineMed_WOOD_FeBC_frac_PD.fillna(0.0)
#FineMed_WOOD_FeBC_frac_PD = FineMed_WOOD_FeBC_frac_PD.where(np.isfinite(FineMed_WOOD_FeBC_frac_PD), 0.0)

#CoaMed_FWOOD_SSP370 = CoaMed_WOOD_FeBC_frac_PD * CoaMed_BCWOOD_SSP370
#CoaMed_FWOOD_SSP370 = CoaMed_FWOOD_SSP370.fillna(0.0)
#CoaMed_FWOOD_SSP370 = CoaMed_FWOOD_SSP370.where(np.isfinite(CoaMed_FWOOD_SSP370), 0.0)

#FineMed_FWOOD_SSP370 = FineMed_WOOD_FeBC_frac_PD * FineMed_BCWOOD_SSP370
#FineMed_FWOOD_SSP370 = FineMed_FWOOD_SSP370.fillna(0.0)
#FineMed_FWOOD_SSP370 = FineMed_FWOOD_SSP370.where(np.isfinite(FineMed_FWOOD_SSP370), 0.0)

#CoaMed_WOOD_FeBC_frac_SSP370 = (CoaMed_FWOOD_SSP370/CoaMed_BCWOOD_SSP370)
#CoaMed_WOOD_FeBC_frac_SSP370 = CoaMed_WOOD_FeBC_frac_SSP370.fillna(0.0)
#CoaMed_WOOD_FeBC_frac_SSP370 = CoaMed_WOOD_FeBC_frac_SSP370.where(np.isfinite(CoaMed_WOOD_FeBC_frac_SSP370), 0.0)

#FineMed_WOOD_FeBC_frac_SSP370 = (FineMed_FWOOD_SSP370/FineMed_BCWOOD_SSP370)
#FineMed_WOOD_FeBC_frac_SSP370 = FineMed_WOOD_FeBC_frac_SSP370.fillna(0.0)
#FineMed_WOOD_FeBC_frac_SSP370 = FineMed_WOOD_FeBC_frac_SSP370.where(np.isfinite(FineMed_WOOD_FeBC_frac_SSP370), 0.0)

# showing ratios of Fe to BC look the same between PD and SSP370 when plotted 
#CoaMed_OIL_FeBC_frac_PD.plot(subplot_kws={'projection': ccrs.PlateCarree()},transform=ccrs.PlateCarree(), cbar_kwargs={'label': 'Ratio'}, vmin =0, vmax =1)
#plt.title("Ratio of Fe to BC PD")
#plt.gca().coastlines()
#plt.show()

#CoaMed_OIL_FeBC_frac_SSP370.plot(subplot_kws={'projection': ccrs.PlateCarree()},transform=ccrs.PlateCarree(), cbar_kwargs={'label': 'Ratio'}, vmin =0, vmax =1)
#plt.title("Ratio of Fe to BC SSP370")
#plt.gca().coastlines()
#plt.show() 

# SMELT ------------------------------------------------------------------------------
CoaMed_FSMELT_SSP370 = CoaMed_FSMELT_PD
CoaMed_FSMELT_SSP370 = CoaMed_FSMELT_SSP370.fillna(0.0)
CoaMed_FSMELT_SSP370 = CoaMed_FSMELT_SSP370.where(np.isfinite(CoaMed_FSMELT_SSP370), 0.0)

FineMed_FSMELT_SSP370 = FineMed_FSMELT_PD
FineMed_FSMELT_SSP370 = FineMed_FSMELT_SSP370.fillna(0.0)
FineMed_FSMELT_SSP370 = FineMed_FSMELT_SSP370.where(np.isfinite(FineMed_FSMELT_SSP370), 0.0)

#CoaMed_FCOALBF_SSP370.plot(subplot_kws={'projection': ccrs.PlateCarree()},transform=ccrs.PlateCarree(), cbar_kwargs={'label': 'Emissions'}, vmin = 0, vmax = 3e-15)
#plt.title("Before zeroes")
#plt.gca().coastlines()
#plt.show() # plot looks right for each sector

# Zeroing out emissions at very end 
CoaMed_FCOALBF_SSP370 = CoaMed_FCOALBF_SSP370.fillna(0.0).where(np.isfinite(CoaMed_FCOALBF_SSP370), 0.0)
CoaMed_FCOALFF_SSP370 = CoaMed_FCOALFF_SSP370.fillna(0.0).where(np.isfinite(CoaMed_FCOALFF_SSP370), 0.0)
CoaMed_FOIL_SSP370 = CoaMed_FOIL_SSP370.fillna(0.0).where(np.isfinite(CoaMed_FOIL_SSP370), 0.0)
CoaMed_FWOOD_SSP370 = CoaMed_FWOOD_SSP370.fillna(0.0).where(np.isfinite(CoaMed_FWOOD_SSP370), 0.0)

FineMed_FCOALBF_SSP370 = FineMed_FCOALBF_SSP370.fillna(0.0).where(np.isfinite(FineMed_FCOALBF_SSP370), 0.0)
FineMed_FCOALFF_SSP370 = FineMed_FCOALFF_SSP370.fillna(0.0).where(np.isfinite(FineMed_FCOALFF_SSP370), 0.0)
FineMed_FOIL_SSP370 = FineMed_FOIL_SSP370.fillna(0.0).where(np.isfinite(FineMed_FOIL_SSP370), 0.0)
FineMed_FWOOD_SSP370 = FineMed_FWOOD_SSP370.fillna(0.0).where(np.isfinite(FineMed_FWOOD_SSP370), 0.0)

CoaMed_BCCOALBF_SSP370 = CoaMed_BCCOALBF_SSP370.fillna(0.0).where(np.isfinite(CoaMed_BCCOALBF_SSP370), 0.0)
CoaMed_BCCOALFF_SSP370 = CoaMed_BCCOALFF_SSP370.fillna(0.0).where(np.isfinite(CoaMed_BCCOALFF_SSP370), 0.0)
CoaMed_BCOIL_SSP370 = CoaMed_BCOIL_SSP370.fillna(0.0).where(np.isfinite(CoaMed_BCOIL_SSP370), 0.0)
CoaMed_BCWOOD_SSP370 = CoaMed_BCWOOD_SSP370.fillna(0.0).where(np.isfinite(CoaMed_BCWOOD_SSP370), 0.0)

FineMed_BCCOALBF_SSP370 = FineMed_BCCOALBF_SSP370.fillna(0.0).where(np.isfinite(FineMed_BCCOALBF_SSP370), 0.0)
FineMed_BCCOALFF_SSP370 = FineMed_BCCOALFF_SSP370.fillna(0.0).where(np.isfinite(FineMed_BCCOALFF_SSP370), 0.0)
FineMed_BCOIL_SSP370 = FineMed_BCOIL_SSP370.fillna(0.0).where(np.isfinite(FineMed_BCOIL_SSP370), 0.0)
FineMed_BCWOOD_SSP370 = FineMed_BCWOOD_SSP370.fillna(0.0).where(np.isfinite(FineMed_BCWOOD_SSP370), 0.0)

# testing final Fe emission plots
#CoaMed_FCOALBF_SSP370.plot(subplot_kws={'projection': ccrs.PlateCarree()},transform=ccrs.PlateCarree(), cbar_kwargs={'label': 'Emissions'}, vmin = 0, vmax = 3e-15)
#plt.title("XX")
#plt.gca().coastlines()
#plt.show() # plot looks right for each sector

combined_Fe_emiss_SSP370 = xr.Dataset({
    "CoaMed_FCOALBF": CoaMed_FCOALBF_SSP370,
    "FineMed_FCOALBF": FineMed_FCOALBF_SSP370,
    "CoaMed_FCOALFF": CoaMed_FCOALFF_SSP370,
    "FineMed_FCOALFF": FineMed_FCOALFF_SSP370,
    "CoaMed_FOIL": CoaMed_FOIL_SSP370,
    "FineMed_FOIL": FineMed_FOIL_SSP370,
    "CoaMed_FWOOD": CoaMed_FWOOD_SSP370,
    "FineMed_FWOOD": FineMed_FWOOD_SSP370,
    "CoaMed_FSMELT": CoaMed_FSMELT_SSP370,
    "FineMed_FSMELT": FineMed_FSMELT_SSP370,
})

combined_Fe_emiss_SSP370["CoaMed_FCOALBF"].attrs["long_name"] = "Coarse Mode Fe from Biofuels"
combined_Fe_emiss_SSP370["CoaMed_FCOALBF"].attrs["units"] = "kg m-2 s-1"
combined_Fe_emiss_SSP370["FineMed_FCOALBF"].attrs["long_name"] = "Fine Mode Fe from Biofuels"
combined_Fe_emiss_SSP370["FineMed_FCOALBF"].attrs["units"] = "kg m-2 s-1"

combined_Fe_emiss_SSP370["CoaMed_FCOALFF"].attrs["long_name"] = "Coarse Mode Fe from Fossil fuels"
combined_Fe_emiss_SSP370["CoaMed_FCOALFF"].attrs["units"] = "kg m-2 s-1"
combined_Fe_emiss_SSP370["FineMed_FCOALFF"].attrs["long_name"] = "Fine Mode Fe from Fossil fuels"
combined_Fe_emiss_SSP370["FineMed_FCOALFF"].attrs["units"] = "kg m-2 s-1"

combined_Fe_emiss_SSP370["CoaMed_FOIL"].attrs["long_name"] = "Coarse Mode Fe from Oil"
combined_Fe_emiss_SSP370["CoaMed_FOIL"].attrs["units"] = "kg m-2 s-1"
combined_Fe_emiss_SSP370["FineMed_FOIL"].attrs["long_name"] = "Fine Mode Fe from Oil"
combined_Fe_emiss_SSP370["FineMed_FOIL"].attrs["units"] = "kg m-2 s-1"

combined_Fe_emiss_SSP370["CoaMed_FWOOD"].attrs["long_name"] = "Coarse Mode Fe from Wood"
combined_Fe_emiss_SSP370["CoaMed_FWOOD"].attrs["units"] = "kg m-2 s-1"
combined_Fe_emiss_SSP370["FineMed_FWOOD"].attrs["long_name"] = "Fine Mode Fe from Wood"
combined_Fe_emiss_SSP370["FineMed_FWOOD"].attrs["units"] = "kg m-2 s-1"

combined_Fe_emiss_SSP370["CoaMed_FSMELT"].attrs["long_name"] = "Coarse Mode Fe from Smelting"
combined_Fe_emiss_SSP370["CoaMed_FSMELT"].attrs["units"] = "kg m-2 s-1"
combined_Fe_emiss_SSP370["FineMed_FSMELT"].attrs["long_name"] = "Fine Mode Fe from Smelting"
combined_Fe_emiss_SSP370["FineMed_FSMELT"].attrs["units"] = "kg m-2 s-1"

combined_Fe_emiss_PD = xr.Dataset({
    "CoaMed_FCOALBF": CoaMed_FCOALBF_PD,
    "FineMed_FCOALBF": FineMed_FCOALBF_PD,
    "CoaMed_FCOALFF": CoaMed_FCOALFF_PD,
    "FineMed_FCOALFF": FineMed_FCOALFF_PD,
    "CoaMed_FOIL": CoaMed_FOIL_PD,
    "FineMed_FOIL": FineMed_FOIL_PD,
    "CoaMed_FWOOD": CoaMed_FWOOD_PD,
    "FineMed_FWOOD": FineMed_FWOOD_PD,
    "CoaMed_FSMELT": CoaMed_FSMELT_PD,
    "FineMed_FSMELT": FineMed_FSMELT_PD,
})

combined_BC_emiss_SSP370 = xr.Dataset({
    "CoaMed_BCCOALBF": CoaMed_BCCOALBF_SSP370,
    "FineMed_BCCOALBF": FineMed_BCCOALBF_SSP370,
    "CoaMed_BCCOALFF": CoaMed_BCCOALFF_SSP370,
    "FineMed_BCCOALFF": FineMed_BCCOALFF_SSP370,
    "CoaMed_BCOIL": CoaMed_BCOIL_SSP370,
    "FineMed_BCOIL": FineMed_BCOIL_SSP370,
    "CoaMed_BCWOOD": CoaMed_BCWOOD_SSP370,
    "FineMed_BCWOOD": FineMed_BCWOOD_SSP370
})

combined_BC_emiss_PD = xr.Dataset({
    "CoaMed_BCCOALBF": CoaMed_BCCOALBF_PD,
    "FineMed_BCCOALBF": FineMed_BCCOALBF_PD,
    "CoaMed_BCCOALFF": CoaMed_BCCOALFF_PD,
    "FineMed_BCCOALFF": FineMed_BCCOALFF_PD,
    "CoaMed_BCOIL": CoaMed_BCOIL_PD,
    "FineMed_BCOIL": FineMed_BCOIL_PD,
    "CoaMed_BCWOOD": CoaMed_BCWOOD_PD,
    "FineMed_BCWOOD": FineMed_BCWOOD_PD
})

# evaluating regions where FeFU might exceed BCFU
# test for each variable and also plot where exceeds for easy visualization
#exceedance_indices = []
#xarr1 = combined_Fe_emiss_PD["CoaMed_FOIL"]
#xarr2 = combined_BC_emiss_PD["CoaMed_BCOIL"]
#for i in range(xarr1.shape[0]):
#    for j in range(xarr1.shape[1]):
#        if xarr1[i, j] > xarr2[i, j]:
#            exceedance_indices.append((i, j))
#print("Exceedance Indices:", exceedance_indices)
#total_indices = xarr1.size  # Total number of elements in the xarrays
#exceedance_count = len(exceedance_indices)
#exceedance_percentage = (exceedance_count / total_indices) * 100
#print(f"Percentage of Exceedance: {exceedance_percentage:.2f}%")
#exceedance_map = (xarr1 > xarr2).astype(int)
#lat = xarr1["lat"]
# #lon = xarr1["lon"]
#import cartopy.feature as cfeature
#import matplotlib.colors as mcolors
#custom_cmap = mcolors.ListedColormap(["white", "red"])
#fig, ax = plt.subplots(
#    figsize=(10, 8),
#    subplot_kw={"projection": ccrs.PlateCarree()})
#im = ax.pcolormesh(
#    lon, lat, exceedance_map.values,
#    cmap=custom_cmap,  # Flipped colormap
#    shading="auto",
#    transform=ccrs.PlateCarree())
#ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
#ax.add_feature(cfeature.BORDERS, linestyle=":")
#gl = ax.gridlines(draw_labels=True, linestyle="--", color="gray", alpha=0.5)
#gl.top_labels = False
#gl.right_labels = False
#ax.set_title("grid cells where COARSE OIL Fe exceeds BC (Present Day)", fontsize=14)
#plt.show()

print(combined_Fe_emiss_SSP370)

# reading out new netcdf file with SSP370 Fe emissions
combined_Fe_emiss_SSP370.to_netcdf("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\coalflyash_FeComb_emissions_SSP370_2040-2050.nc")
# This looks good in Panoply, no apparent gridding errors, Now to calculate budgets and compare future to PD emissions

# COMPARING BC vs Fe in SSP370 -- BC does look higher 
CoaMed_BCOIL_SSP370.plot(subplot_kws={'projection': ccrs.PlateCarree()},transform=ccrs.PlateCarree(), 
cbar_kwargs={'label': 'emission'}, vmin=1e-18,vmax=1e-15)
plt.title("Black Carbon Fu")
plt.gca().coastlines()
plt.show() 

CoaMed_FOIL_SSP370.plot(subplot_kws={'projection': ccrs.PlateCarree()},transform=ccrs.PlateCarree(), cbar_kwargs={'label': 'emission'}, vmin=1e-18,vmax=1e-15)
plt.title("Iron Fu")
plt.gca().coastlines()
plt.show() 

CoaMed_BCOIL_PD.plot(subplot_kws={'projection': ccrs.PlateCarree()},transform=ccrs.PlateCarree(), 
cbar_kwargs={'label': 'emission'}, vmin=1e-18,vmax=1e-15)
plt.title("Black Carbon PD")
plt.gca().coastlines()
plt.show() 

CoaMed_FOIL_PD.plot(subplot_kws={'projection': ccrs.PlateCarree()},transform=ccrs.PlateCarree(), cbar_kwargs={'label': 'emission'}, vmin=1e-18,vmax=1e-15)
plt.title("Iron PD")
plt.gca().coastlines()
plt.show() 
# oil emissions lower in SSP370 scenario ? 

#%% CALCULATING EMISSION BUDGETS -------------------------------
import numpy as np
import xarray as xr

# re-adding iron together for final budget ratio checks 
all_variables = xr.Dataset({
    "CoaMed_FCOALBF_SSP370": CoaMed_FCOALBF_SSP370,
    "FineMed_FCOALBF_SSP370": FineMed_FCOALBF_SSP370,
    "CoaMed_FCOALFF_SSP370": CoaMed_FCOALFF_SSP370,
    "FineMed_FCOALFF_SSP370": FineMed_FCOALFF_SSP370,
    "CoaMed_FOIL_SSP370": CoaMed_FOIL_SSP370,
    "FineMed_FOIL_SSP370": FineMed_FOIL_SSP370,
    "CoaMed_FWOOD_SSP370": CoaMed_FWOOD_SSP370,
    "FineMed_FWOOD_SSP370": FineMed_FWOOD_SSP370,
    "CoaMed_FSMELT_SSP370": CoaMed_FSMELT_SSP370,
    "FineMed_FSMELT_SSP370": FineMed_FSMELT_SSP370,

    "CoaMed_FCOALBF_PD": CoaMed_FCOALBF_PD,
    "FineMed_FCOALBF_PD": FineMed_FCOALBF_PD,
    "CoaMed_FCOALFF_PD": CoaMed_FCOALFF_PD,
    "FineMed_FCOALFF_PD": FineMed_FCOALFF_PD,
    "CoaMed_FOIL_PD": CoaMed_FOIL_PD,
    "FineMed_FOIL_PD": FineMed_FOIL_PD,
    "CoaMed_FWOOD_PD": CoaMed_FWOOD_PD,
    "FineMed_FWOOD_PD": FineMed_FWOOD_PD,
    "CoaMed_FSMELT_PD": CoaMed_FSMELT_PD,
    "FineMed_FSMELT_PD": FineMed_FSMELT_PD,

    "CoaMed_BCCOALBF_SSP370": CoaMed_BCCOALBF_SSP370,
    "FineMed_BCCOALBF_SSP370": FineMed_BCCOALBF_SSP370,
    "CoaMed_BCCOALFF_SSP370": CoaMed_BCCOALFF_SSP370,
    "FineMed_BCCOALFF_SSP370": FineMed_BCCOALFF_SSP370,
    "CoaMed_BCOIL_SSP370": CoaMed_BCOIL_SSP370,
    "FineMed_BCOIL_SSP370": FineMed_BCOIL_SSP370,
    "CoaMed_BCWOOD_SSP370": CoaMed_BCWOOD_SSP370,
    "FineMed_BCWOOD_SSP370": FineMed_BCWOOD_SSP370,

    "CoaMed_BCCOALBF_PD": CoaMed_BCCOALBF_PD,
    "FineMed_BCCOALBF_PD": FineMed_BCCOALBF_PD,
    "CoaMed_BCCOALFF_PD": CoaMed_BCCOALFF_PD,
    "FineMed_BCCOALFF_PD": FineMed_BCCOALFF_PD,
    "CoaMed_BCOIL_PD": CoaMed_BCOIL_PD,
    "FineMed_BCOIL_PD": FineMed_BCOIL_PD,
    "CoaMed_BCWOOD_PD": CoaMed_BCWOOD_PD,
    "FineMed_BCWOOD_PD": FineMed_BCWOOD_PD,

    "BothMode_FCOALBF_SSP370":BothMode_FCOALBF_SSP370,
    "BothMode_FCOALFF_SSP370":BothMode_FCOALFF_SSP370,
    "BothMode_FOIL_SSP370":BothMode_FOIL_SSP370,
    "BothMode_FWOOD_SSP370":BothMode_FWOOD_SSP370,

    "BothMode_FCOALBF_PD":BothMode_FCOALBF_PD,
    "BothMode_FCOALFF_PD":BothMode_FCOALFF_PD,
    "BothMode_FOIL_PD":BothMode_FOIL_PD,
    "BothMode_FWOOD_PD":BothMode_FWOOD_PD,

    "BothMode_BC_COALBF_SSP370":BC_COALBF_SSP370,
    "BothMode_BC_COALFF_SSP370":BC_COALFF_SSP370,
    "BothMode_BC_OIL_SSP370":BC_OIL_SSP370,
    "BothMode_BC_WOOD_SSP370":BC_WOOD_SSP370,

    "BothMode_BC_COALBF_PD":BC_COALBF_PD,
    "BothMode_BC_COALFF_PD":BC_COALFF_PD,
    "BothMode_BC_OIL_PD":BC_OIL_PD,
    "BothMode_BC_WOOD_PD":BC_WOOD_PD,
})

ds_file = all_variables

#ds_file = combined_Fe_emiss_SSP370

#ds_file = combined_Fe_emiss_PD

#ds_file = combined_BC_emiss_SSP370

#ds_file = combined_BC_emiss_PD

#print(ds_file)

# online way to calculate emission budgets 
# Load latitude and longitude values
lon = ds_file['lon'].values  # Longitude in degrees
lat = ds_file['lat'].values  # Latitude in degrees

# Convert to radians
lon_rads = np.radians(lon)
lat_rads = np.radians(lat)

# Compute grid spacings (assuming equidistant grid)
d_lat = np.abs(lat[1] - lat[0])  # Latitude grid spacing in degrees
d_lon = np.abs(lon[1] - lon[0])  # Longitude grid spacing in degrees

# Half-grid spacing for staggered grid
g_lat = np.radians(d_lat / 2)  # Latitude half-spacing in radians
g_lon = np.radians(d_lon / 2)  # Longitude half-spacing in radians

# Radius of Earth in m
R = 6.3781E6 

# Empty list to store results
cell_areas_staggered = []

# Loop through each grid cell
# specifies the coordinates with [j,i] -- lat being i and lon being j 
for i in range(len(lat)):
    for j in range(len(lon)):
        # Convert latitude and longitude to radians
        lat_center = lat_rads[i]
        lon_center = lon_rads[j]

        # Compute staggered latitudes per Arakawa Lamb - C gridding for CESM
        lat_north = lat_center + g_lat
        lat_south = lat_center - g_lat

        # Ensure staggered latitudes are within valid range (-π/2 to π/2)
        lat_north = np.clip(lat_north, -np.pi / 2, np.pi / 2)
        lat_south = np.clip(lat_south, -np.pi / 2, np.pi / 2)

        # Compute area of the cell
        area = R**2 * (np.sin(lat_north) - np.sin(lat_south)) * (2 * g_lon)
        cell_areas_staggered.append(area)

# Convert to numpy array and reshape to match grid dimensions
cell_areas_staggered = np.array(cell_areas_staggered).reshape(len(lat), len(lon))

# Verify to see if areas add to 5.1E14 in m^2
sum_sa_earth = cell_areas_staggered.sum()
print(f"surface area, {sum_sa_earth:3e}") 

# add cell area to original xarray for easy calling 
ds_file['cell_area'] = xr.DataArray(
    cell_areas_staggered,
    dims=["lat", "lon"],  # Same dimensions as in the original dataset
    coords={"lat": ds_file['lat'], "lon": ds_file['lon']},  # Use original coordinates
    attrs={
        "units": "m^2",  # Specify units for the cell area
        "description": "Calculated grid cell area using staggered grid approach",
    },
)

print(ds_file)

# Calculating budgets-- started with OG Fe emissions to check against MIMI, got the same results as the for loop 
import pandas as pd
# Dictionaries to store both individual and total budgets
global_budgets = {}
individual_cell_emissions = {}
# Not sure if this loop is running properly, getting different results when I do it manually
# Loop over all variables in the dataset
for var_name in ds_file.data_vars:
    if "Mode" in var_name:  
        # Calculate individual emissions (for each cell) and convert from sec to annual and from kg to Tg
        # Assuming input variable is the flux in kg/m²/s and 'cell_area' is the area of each grid cell in m²
        individual_cell_emissions[var_name] = (ds_file[var_name] * ds_file['cell_area'] * 3600 * 24 * 365 *1E-9)
        
        # Calculate the total budget by summing individual emissions
        total_budget = individual_cell_emissions[var_name].sum()
        
        # Store the total budget in the dictionary
        global_budgets[var_name] = total_budget
        
        # Print the budget for the variable
        print(f"FU Total Fe budget of {var_name} (Tg/year): {total_budget:.3e}")

# Create a DataFrame from the budgets dictionary (Total Budgets)
budget_df = pd.DataFrame(list(global_budgets.items()), columns=["Variable", "Total_Budget_Tg_per_year"])
print(budget_df)

## Create a DataFrame for the individual emissions, flattening the xarray values
#individual_df = pd.DataFrame()
#for var_name in individual_cell_emissions:
    # Flatten the DataArray into a pandas DataFrame for the individual emissions
    #individual_df[var_name] = individual_cell_emissions[var_name].values.flatten()

# Save the total budgets and individual emissions to separate sheets in an Excel file
output_file = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\Budgets_01282025.xlsx"
with pd.ExcelWriter(output_file) as writer:
   budget_df.to_excel(writer, sheet_name="Total_Budgets", index=False)  
   print(f"Total budgets and individual emissions saved to {output_file}")

#%%
ds_file['CoaMed_FOIL_SSP370_Tg'] = ds_file['CoaMed_FOIL_SSP370']*ds_file['cell_area']*3600*24*365*1E-9
ds_file['CoaMed_FOIL_PD_Tg'] = ds_file['CoaMed_FOIL_PD']*ds_file['cell_area']*3600*24*365*1E-9
ds_file['CoaMed_BCOIL_SSP370_Tg'] = ds_file['CoaMed_BCOIL_SSP370']*ds_file['cell_area']*3600*24*365*1E-9
ds_file['CoaMed_BCOIL_PD_Tg'] = ds_file['CoaMed_BCOIL_PD']*ds_file['cell_area']*3600*24*365*1E-9
FePD_FeFU = ds_file['CoaMed_FOIL_PD_Tg']/ds_file['CoaMed_FOIL_SSP370_Tg']
FePD_FeFU = FePD_FeFU.fillna(0.0) # this didnt seem to make a difference 
FePD_FeFU = FePD_FeFU.where(np.isfinite(FePD_FeFU), 0.0)
ds_file['FePD:FeFU'] = FePD_FeFU
BCPD_BCFU = ds_file['CoaMed_BCOIL_PD_Tg']/ds_file['CoaMed_BCOIL_SSP370_Tg']
BCPD_BCFU = BCPD_BCFU.fillna(0.0)
BCPD_BCFU = BCPD_BCFU.where(np.isfinite(BCPD_BCFU), 0.0)
ds_file['BCPD:BCFU'] = BCPD_BCFU

CoaMed_FOIL_SSP370_Tg = ds_file['CoaMed_FOIL_SSP370']*ds_file['cell_area']*3600*24*365*1E-9
total_Fe_SSP370_emiss = CoaMed_FOIL_SSP370_Tg.sum()
print("Fe_SSP370_Tg", total_Fe_SSP370_emiss)

CoaMed_FOIL_PD_Tg = ds_file['CoaMed_FOIL_PD']*ds_file['cell_area']*3600*24*365*1E-9
total_Fe_PD_emiss = CoaMed_FOIL_PD_Tg.sum()
print("Fe_PD_Tg", total_Fe_PD_emiss)

CoaMed_BCOIL_SSP370_Tg = ds_file['CoaMed_BCOIL_SSP370']*ds_file['cell_area']*3600*24*365*1E-9
total_BC_SSP370_emiss = CoaMed_BCOIL_SSP370_Tg.sum()
print("BC_SSP370_Tg", total_BC_SSP370_emiss)

CoaMed_BCOIL_PD_Tg = ds_file['CoaMed_BCOIL_PD']*ds_file['cell_area']*3600*24*365*1E-9
total_BC_PD_emiss = CoaMed_BCOIL_PD_Tg.sum()
print("BC_PD_Tg", total_BC_PD_emiss)

quick_total_ratio_PD = total_Fe_PD_emiss/total_BC_PD_emiss
print(quick_total_ratio_PD)
quick_total_ratio_SSP370 = total_Fe_SSP370_emiss/total_BC_SSP370_emiss
print(quick_total_ratio_SSP370)

#print(ds_file)

ds_file['CoaMed_FWOOD_SSP370_Tg'] = ds_file['CoaMed_FWOOD_SSP370']*ds_file['cell_area']*3600*24*365*1E-9
ds_file['CoaMed_FWOOD_PD_Tg'] = ds_file['CoaMed_FWOOD_PD']*ds_file['cell_area']*3600*24*365*1E-9
ds_file['CoaMed_BCWOOD_SSP370_Tg'] = ds_file['CoaMed_BCWOOD_SSP370']*ds_file['cell_area']*3600*24*365*1E-9
ds_file['CoaMed_BCWOOD_PD_Tg'] = ds_file['CoaMed_BCWOOD_PD']*ds_file['cell_area']*3600*24*365*1E-9
FePD_FeFU = ds_file['CoaMed_FWOOD_PD_Tg']/ds_file['CoaMed_FWOOD_SSP370_Tg']
FePD_FeFU = FePD_FeFU.fillna(0.0)
FePD_FeFU = FePD_FeFU.where(np.isfinite(FePD_FeFU), 0.0)
ds_file['FePD:FeFU'] = FePD_FeFU
BCPD_BCFU = ds_file['CoaMed_BCWOOD_PD_Tg']/ds_file['CoaMed_BCWOOD_SSP370_Tg']
BCPD_BCFU = BCPD_BCFU.fillna(0.0)
BCPD_BCFU = BCPD_BCFU.where(np.isfinite(BCPD_BCFU), 0.0)
ds_file['BCPD:BCFU'] = BCPD_BCFU

CoaMed_FWOOD_SSP370_Tg = ds_file['CoaMed_FWOOD_SSP370']*ds_file['cell_area']*3600*24*365*1E-9
total_Fe_SSP370_emiss = CoaMed_FWOOD_SSP370_Tg.sum()
print("Fe_SSP370_Tg", total_Fe_SSP370_emiss)

CoaMed_FWOOD_PD_Tg = ds_file['CoaMed_FWOOD_PD']*ds_file['cell_area']*3600*24*365*1E-9
total_Fe_PD_emiss = CoaMed_FWOOD_PD_Tg.sum()
print("Fe_PD_Tg", total_Fe_PD_emiss)

CoaMed_BCWOOD_SSP370_Tg = ds_file['CoaMed_BCWOOD_SSP370']*ds_file['cell_area']*3600*24*365*1E-9
total_BC_SSP370_emiss = CoaMed_BCWOOD_SSP370_Tg.sum()
print("BC_SSP370_Tg", total_BC_SSP370_emiss)

CoaMed_BCWOOD_PD_Tg = ds_file['CoaMed_BCWOOD_PD']*ds_file['cell_area']*3600*24*365*1E-9
total_BC_PD_emiss = CoaMed_BCWOOD_PD_Tg.sum()
print("BC_PD_Tg", total_BC_PD_emiss)

quick_total_ratio_PD = total_Fe_PD_emiss/total_BC_PD_emiss
print(quick_total_ratio_PD)
quick_total_ratio_SSP370 = total_Fe_SSP370_emiss/total_BC_SSP370_emiss
print(quick_total_ratio_SSP370)

# for individual sector calculations -----------------------------------------------------------------------
# conversion from kg/m2/s to kg per year 
# [FOILi,j]*cell_areai,j*3600*24*365 = [FOILi,j]
# sum[FOIL] for all i,j combinations 
# want to store value in kg for each grid cell as mass_FOIL within the netcdf as a new variable but also sum as a budget_FOIL as a single value I can then print

# FUTURE INPUT DATASET (SSP370)
# Assuming 'FOIL' is the flux in kg/m²/s and 'cell_area' is the area of each grid cell in m²
# the .keys() quickly shows the names of the variables in your netcdf
# ds_FU_file.keys()

# Convert FOIL to kg per year for each grid cell (kg)
#ds_FU_file['mass_CoaMed_FOIL'] = ds_FU_file['CoaMed_FOIL'] * ds_FU_file['cell_area'] * 3600 * 24 * 365

# Sum the mass_FOIL over all grid cells to get the total budget in kg
#budget_FOIL = ds_FU_file['mass_CoaMed_FOIL'].sum()
#print(f"Total budget of CoaMed_FOIL (kg/year): {budget_FOIL.values:.3e}")

# Show the DataFrame interactively
#df = ds_FU_file['mass_CoaMed_FOIL'].to_dataframe()
#from itables import show
#df_subset = df.iloc[20000:21000]
#show(df_subset, scrollY=True, scrollX=False, maxRows=100, maxColumns=10)
# confirmed there is data in the cells ----------------------------------------------------------------------------------

# using cdo method to calculate grid areas offline ----------------------------------------------------------------------------------
#grid_area_file = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\192x288_gridarea.nc"
#grid_areas = xr.open_dataset(grid_area_file)
# used cdo grid area to do this
#grid_areas['lon'] = ds_FU_file['lon'] # ensuring lats and lons are identical
#grid_areas['lat'] = ds_FU_file['lat']
#ds_FU_file['cell_area'] = grid_areas['cell_area']
# #gridarea = ds_FU_file['cell_area']  #grid cell areas
#sum_gridarea = gridarea.sum()
#print(sum_gridarea) # 5.1 x 10^14 square meters, passed check of surface area of the earth
# there are subtle differences between the grid areas calculated using cdo gridarea and this script but I think small enough that it is okay -- the poles are the exact same which was my previous issue
# cdo remapbil,grid_dims.txt inputfile.nc outputfile.nc OR cdo -f nc -setgrid,grid_dims.txt inputfile.nc outputfile.nc
# remapcon better for emissions, see differences between bilineation and conservative methods
# ------------------------------------------------------------------------------------------------------------------------

#%% OLD CODE BELOW ---- DONT RUN THIS REGRIDDER
#from netCDF4 import Dataset

# Quick view of some of the details within the netcdf file
# Open the dataset 
#dataset = Dataset('C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\BC-em-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_2000-2015MEAN.nc', 'r')

# List all variables and attributes
#print(dataset.variables.keys())
#print(dataset.variables['BC_em_anthro'])
#print(dataset.variables['sector'])
#print(dataset.variables['sector_bnds'])

# 0: Agriculture; 1: Energy;  2: Industrial; 3: Transportation; 4: Residential, Commercial, Other; 5: Solvents production and application; 6: Waste; 7: International Shipping

# FCOAL = 1 + 2 (FF) + 4 (BF)
# FOIL = 3 (only from terrestrial cells) + 7 (only from ocean cells)
# FWOOD = 4 (BF?)
# FSMELT = bring over from Fe_Emissions_Fuel
# I think I actually need to split coal into the following though, based on model variable pointers
# FCOALFF = 1 + 2
# FCOALBF = 4

# F370[i,j,s] = BC370[i,j,s] / BCPD[i,j,s] * FePD[i,j,s]
# where i is lon, j is lat, and s is sector 
# calculate for each sector, and then add to create new consolidated variables
# need to double check with Douglas on wood and FCOAL

# visualizing the BC_emiss data for each sector in separate pandas dataframes
#import xarray as xr
#import pandas as pd
#import numpy as np

# Load the NetCDF file -- for Future simulation
#file_path = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\BC-em-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_2000-2015MEAN.nc"
#ds = xr.open_dataset(file_path)

# grid dimensions are currently 360 lat x 720 lon
# re-gridding to 192 lat x 288 lon
#import xesmf as xe

# Create target grid (192 x 288)
#lat_target = xr.DataArray(
    #data=np.linspace(-90, 90, 192),
   # dims='lat',
    #name='lat')
#lon_target = xr.DataArray(
  #  data=np.linspace(0, 360, 288, endpoint=False),  # 0 to 360 for 288 points
   # dims='lon',
    #name='lon')
#target_grid = xr.Dataset({'lat': lat_target, 'lon': lon_target})

# Apply regridding
#regridder = xe.Regridder(ds, target_grid, method='bilinear')
#ds_regridded = regridder(ds)

# Save to a new NetCDF file
#ds_regridded.to_netcdf('C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\BC-em-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_2000-2015MEAN_regridded_288x192.nc') 

# now for SSP370 File
#from netCDF4 import Dataset
# Open the dataset 
#dataset = Dataset('C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\BC-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-AIM-ssp370-1-1_gn_2090-2100MEAN.nc', 'r')

# List all variables and attributes
#print(dataset.variables.keys())
#print(dataset.variables['BC_em_anthro'])
#print(dataset.variables['sector'])
#print(dataset.variables['sector_bnds'])

# the time variables needs to be condensed into a mean first -- time is 120 options, and need to regrid
#import xarray as xr

# Load the NetCDF file -- for Future simulation
#file_path = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\BC-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-AIM-ssp370-1-1_gn_2090_2100MEAN.nc"
#ds = xr.open_dataset(file_path)
#BC_em_anthro = ds['BC_em_anthro']

# Calculate the mean across the `time` dimension
#mean_BC_em_anthro = BC_em_anthro.mean(dim="time", keep_attrs=True)

#mean_BC_em_anthro.to_netcdf("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\BC-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-AIM-ssp370-1-1_gn_MEAN.nc")

# ----------------------------------------------------------------------------------------------------------------
# grid dimensions are currently 360 lat x 720 lon
# re-gridding to 192 lat x 288 lon
#import xarray as xr
#import xesmf as xe
#import numpy as np

#file_path = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\BC-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-AIM-ssp370-1-1_gn_2090_2100MEAN.nc"
#ds = xr.open_dataset(file_path)


# Create target grid (192 x 288)
#lat_target = xr.DataArray(
   # data=np.linspace(-90, 90, 192),
   # dims='lat',
    #name='lat')
#lon_target = xr.DataArray(
   # data=np.linspace(0, 360, 288, endpoint=False),  # 0 to 360 for 288 points
   ## dims='lon',
   # name='lon')
#target_grid = xr.Dataset({'lat': lat_target, 'lon': lon_target})

# Apply regridding
#regridder = xe.Regridder(ds, target_grid, method='bilinear')
#ds_regridded = regridder(ds)

# Save to a new NetCDF file
#ds_regridded.to_netcdf('C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\BC-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-AIM-ssp370-1-1_gn_2090-2100MEAN_regridded_288x192.nc')

# this whole regridder didnt work, 
# instead it is necessary to use cdo remapbil 

# %% DATAFRAMES FOR QUICK NUMBERS VIZ IN TABLE --------------------------------------------------------------------------
# visualizing the BC_emiss data for each sector in separate pandas dataframes
import xarray as xr
import pandas as pd
import numpy as np

# Load the NetCDF file -- using newly regridded data
file_path = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\BC-em-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_2000-2015MEAN_regridded_288x192.nc"
ds = xr.open_dataset(file_path)

# Extract Black Carbon (BC) emission values SSP370 
BC_em_anthro = ds['BC_em_anthro']  # Emission data
lon = ds['lon'].values           # Longitude values -- the .values extracts as a numPy array automatically, removes metadata
lat = ds['lat'].values            # Latitude values
sector = ds['sector'].values      # Sector identifiers

# Check
print("lon", lon)
print("lat", lat)
print("sector", sector) # looks right 

# dataframes are nice for quick viz but not great for manipulation and plotting, xarray is best 
# Create a dictionary to store dataframes for each sector
sector_dataframes = {}

# For loop to loop through individual sectors
for sec_idx, sec_value in enumerate(sector):
    # Filter data for the current sector index
    filtered_data = BC_em_anthro.isel(sector=sec_idx).squeeze()
    
    # Convert the data to a pandas DataFrame
    df = pd.DataFrame(filtered_data.values, index=lat, columns=lon)
    df.index.name = 'Latitude'
    df.columns.name = 'Longitude'
    
    # Store the DataFrame in the dictionary
    sector_dataframes[sec_value] = df

# Access the dataframe for a specific sector, e.g., sector at midpoint `sector_value`
BC_emiss_0 = sector[0]  # Adjust as needed
BC_emiss_0_df = sector_dataframes[BC_emiss_0]
# BC_emiss_0_df # pop out to jupyter notebook 

# Sector 1: Energy
BC_emiss_1 = sector[1]  
BC_emiss_1_df = sector_dataframes[BC_emiss_1]
# Sector 2: Industrial Coal
BC_emiss_2 = sector[2] 
BC_emiss_2_df = sector_dataframes[BC_emiss_2]
# Sector 3: Transportation
BC_emiss_3 = sector[3]  
BC_emiss_3_df = sector_dataframes[BC_emiss_3]
# Sector 4: Residential Coal
BC_emiss_4 = sector[4]  # Adjust as needed
BC_emiss_4_df = sector_dataframes[BC_emiss_4]
# Sector 5: Solvents, etc. 
BC_emiss_5 = sector[5]  # Adjust as needed
BC_emiss_5_df = sector_dataframes[BC_emiss_5]
# Sector 6: Waste
BC_emiss_6 = sector[6]  # Adjust as needed
BC_emiss_6_df = sector_dataframes[BC_emiss_6]
# Sector 7: International Shipping
BC_emiss_7 = sector[7]  # Adjust as needed
BC_emiss_7_df = sector_dataframes[BC_emiss_7]

# used this to do a quick check to make sure different values were coming up uniquely for each sector
from itables import show
# Show the DataFrame interactively
show(BC_emiss_0_df, scrollY=True, scrollX=True, maxRows=10, maxColumns=200) # looks good, data is dif between sectors
print(lon.shape, "lon dim")
print(lat.shape, "lat values")






