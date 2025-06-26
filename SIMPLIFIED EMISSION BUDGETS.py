#%% SIMPLIFIED EMISSION BUDGETS -- For regional and global budgets 
import numpy as np
import xarray as xr
from netCDF4 import Dataset

# LOAD IN THE NETCDF FILE -----------------------------------------
# PI FILE STRUCTURE
#sim_directory = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\CAM6-MIMI-SSP370-INDCOAL0.2-RESICOAL33-WOOD56-OIL38-FIRE33.cam.h1.2009-2011.nc"
# PD FILE STRUCTURE
#sim_directory ="C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\CAM6-MIMI_2010CLIMO_INDCOAL0.2-RESICOAL33-WOOD56-OIL38_2009-2011.h1.BBx1_ANx2_DUx2.nc"

sim_directory = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\CAM6-MIMI_2010CLIMO_INDCOAL0.05-RESICOAL33-WOOD56-OIL25.cam.h1.2009-2011_PD.nc"
ds_h1 = xr.open_dataset(sim_directory)

# EXTRACT VARIABLES ---------------------------------------------------------------
FESOLDRY_var = ds_h1["FESOLDRY"]
FESOLWET_var = ds_h1["FESOLWET"]
FESOLDEP = FESOLDRY_var + FESOLWET_var
FESOLDEP_mean = FESOLDEP.mean(dim="time") # turn this on when data is not yet aggregated over time
#FESOLDEP_mean = FESOLDEP_mean.fillna(0.0).where(np.isfinite(FESOLDEP_mean), 0.0)

FEANSOLDRY_var = ds_h1["FEANSOLDRY"]
FEANSOLWET_var = ds_h1["FEANSOLWET"]
FEANSOLDEP = FEANSOLDRY_var + FEANSOLWET_var
FEANSOLDEP_mean = FEANSOLDEP.mean(dim="time")
#FEANSOLDEP_mean = FEANSOLDEP_mean.fillna(0.0).where(np.isfinite(FEANSOLDEP_mean), 0.0)

FEBBSOLDRY_var = ds_h1["FEBBSOLDRY"]
FEBBSOLWET_var = ds_h1["FEBBSOLWET"]
FEBBSOLDEP = FEBBSOLDRY_var + FEBBSOLWET_var
FEBBSOLDEP_mean = FEBBSOLDEP.mean(dim="time")
#FEBBSOLDEP_mean = FEBBSOLDEP_mean.fillna(0.0).where(np.isfinite(FEBBSOLDEP_mean), 0.0)

FEDUSOLDRY_var = ds_h1["FEDUSOLDRY"]
FEDUSOLWET_var = ds_h1["FEDUSOLWET"]
FEDUSOLDEP = FEDUSOLDRY_var + FEDUSOLWET_var
FEDUSOLDEP_mean = FEDUSOLDEP.mean(dim="time")
#FEDUSOLDEP_mean = FEDUSOLDEP_mean.fillna(0.0).where(np.isfinite(FEDUSOLDEP_mean), 0.0)

FETOTDRY_var = ds_h1["FETOTDRY"]
FETOTWET_var = ds_h1["FETOTWET"]
FETOTDEP = FETOTDRY_var + FETOTWET_var
FETOTDEP_mean = FETOTDEP.mean(dim="time")
#FETOTDEP_mean = FETOTDEP_mean.fillna(0.0).where(np.isfinite(FETOTDEP_mean), 0.0)

FEANTOTDRY_var = ds_h1["FEANTOTDRY"]
FEANTOTWET_var = ds_h1["FEANTOTWET"]
FEANTOTDEP = FEANTOTDRY_var + FEANTOTWET_var
FEANTOTDEP_mean = FEANTOTDEP.mean(dim="time")
#FEANTOTDEP_mean = FEANTOTDEP_mean.fillna(0.0).where(np.isfinite(FEANTOTDEP_mean), 0.0)

FEBBTOTDRY_var = ds_h1["FEBBTOTDRY"]
FEBBTOTWET_var = ds_h1["FEBBTOTWET"]
FEBBTOTDEP = FEBBTOTDRY_var + FEBBTOTWET_var
FEBBTOTDEP_mean = FEBBTOTDEP.mean(dim="time")
#FEBBTOTDEP_mean = FEBBTOTDEP_mean.fillna(0.0).where(np.isfinite(FEBBTOTDEP_mean), 0.0)

FEDUTOTDRY_var = ds_h1["FEDUTOTDRY"]
FEDUTOTWET_var = ds_h1["FEDUTOTWET"]
FEDUTOTDEP = FEDUTOTDRY_var + FEDUTOTWET_var
FEDUTOTDEP_mean = FEDUTOTDEP.mean(dim="time")
#FEDUTOTDEP_mean = FEDUTOTDEP_mean.fillna(0.0).where(np.isfinite(FEDUTOTDEP_mean), 0.0)

# ADD ALL VARIABLES TO SINGLE DATASET --------------------------------------------
all_FE_DEP = xr.Dataset({
    # TOTAL SOURCES
    "FESOLDEP_mean": FESOLDEP_mean,
    "FEANSOLDEP_mean": FEANSOLDEP_mean,
    "FEBBSOLDEP_mean": FEBBSOLDEP_mean,
    "FEDUSOLDEP_mean": FEDUSOLDEP_mean,
    
    "FETOTDEP_mean": FETOTDEP_mean,
    "FEANTOTDEP_mean": FEANTOTDEP_mean,
    "FEBBTOTDEP_mean": FEBBTOTDEP_mean,
    "FEDUTOTDEP_mean": FEDUTOTDEP_mean,
    })


# FIND AREA OF INDIVIDUAL GRID CELLS --------------------------------
lon = all_FE_DEP['lon'].values  # Longitude in degrees
lat = all_FE_DEP['lat'].values  # Latitude in degrees

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

# add cell area to original xarrays for easy calling 
all_FE_DEP['cell_area'] = xr.DataArray(
    cell_areas_staggered,
    dims=["lat", "lon"],  # Same dimensions as in the original dataset
    coords={"lat": all_FE_DEP['lat'], "lon": all_FE_DEP['lon']},  # Use original coordinates
    attrs={
        "units": "m^2",  # Specify units for the cell area
        "description": "Calculated grid cell area using staggered grid approach",
    },
)

# ---- ASSIGNING ATM. DEPOSITION REGIONS 
import pandas as pd
# Dictionaries to store both individual and total budgets
global_budgets = {}
individual_cell_emissions = {}

# Define regional conditions
# DO NOT CHANGE THESE!!! THEY ARE GOOD NOW FOR ADDING FINAL BUDGETS TO TOTAL 
# South Atlantic
SATL_condition = (all_FE_DEP['lon'] >= 295.0) & (all_FE_DEP['lat'] <= -15.0) & (all_FE_DEP['lat'] > -48.0) | \
                 (all_FE_DEP['lon'] < 30) & (all_FE_DEP['lat'] <= -15.0) & (all_FE_DEP['lat'] >-48.0)
# North Atlantic
NATL_condition = (all_FE_DEP['lon'] > 265.0) & (all_FE_DEP['lat'] > 29.0) & (all_FE_DEP['lat'] < 60.0) | \
                 (all_FE_DEP['lon'] < 110.0) & (all_FE_DEP['lat'] > 29.0) & (all_FE_DEP['lat'] < 60.0)
# Arabian Sea
AS_condition = (all_FE_DEP['lon'] < 78.0) & (all_FE_DEP['lon'] >= 30.0) & \
               (all_FE_DEP['lat'] <= 29.0) & (all_FE_DEP['lat'] > -48.0)
# Bay of Bengal
BB_condition = (all_FE_DEP['lon'] < 110.0) & (all_FE_DEP['lon'] >= 78.0) & \
               (all_FE_DEP['lat'] <= 29.0) & (all_FE_DEP['lat'] > -48.0)
# Indian Ocean 
INDO_condition = (all_FE_DEP['lon'] < 110.0) & (all_FE_DEP['lon'] >= 30.0) & \
               (all_FE_DEP['lat'] <= 29.0) & (all_FE_DEP['lat'] > -48.0)
# Southeastern Asia
SEAS_condition = (all_FE_DEP['lon'] < 150.0) & (all_FE_DEP['lon'] >= 110.0) & \
               (all_FE_DEP['lat'] < 60.0) & (all_FE_DEP['lat'] > -15.0)
# ENorth Pacific
ENPAC_condition = (all_FE_DEP['lon'] <= 180.0) & (all_FE_DEP['lon'] >= 150.0) & \
               (all_FE_DEP['lat'] < 60.0) & (all_FE_DEP['lat'] > 29.0)
# WNorth Pacific
WNPAC_condition = (all_FE_DEP['lon'] <= 265.0) & (all_FE_DEP['lon'] >= 180.0) & \
               (all_FE_DEP['lat'] < 60.0) & (all_FE_DEP['lat'] > 29.0)
# North Pacific
NPAC_condition = (all_FE_DEP['lon'] <= 265.0) & (all_FE_DEP['lon'] >= 150.0) & \
               (all_FE_DEP['lat'] < 60.0) & (all_FE_DEP['lat'] > 29.0)
# Arctic
ARCT_condition = (all_FE_DEP['lon'] >= 0.0) & (all_FE_DEP['lat'] >= 60.0) &\
               (all_FE_DEP['lon'] >= 0.0)  & (all_FE_DEP['lat'] <= 90.0)
# Australia/South Pacific
AUSP_condition = (all_FE_DEP['lon'] < 295.0) & (all_FE_DEP['lon'] >= 110.0) & \
               (all_FE_DEP['lat'] <= -15.0) & (all_FE_DEP['lat'] > -48.0)
# Southern Ocean
SO_condition = (all_FE_DEP['lon'] >= 0.0) & (all_FE_DEP['lat'] <= -48.0) &\
               (all_FE_DEP['lon'] >= 0.0)  & (all_FE_DEP['lat'] >= -90.0)
# Central Pacific / Asia
CPAO_condition = (all_FE_DEP['lon'] >= 150.0) & (all_FE_DEP['lat'] > -15.0) & (all_FE_DEP['lat'] <= 29.0) | \
                 (all_FE_DEP['lon'] < 30.0) & (all_FE_DEP['lat'] > -15.0) & (all_FE_DEP['lat'] <= 29.0)

# Initialize an empty dictionary to store the grid cell counts
region_cell_counts = {}
total_cell_counts = {}

# Loop through each region condition and count the grid cells that satisfy the condition
region_cell_counts['South Atlantic'] = (SATL_condition).sum().values
region_cell_counts['North Atlantic'] = (NATL_condition).sum().values
region_cell_counts['Arabian Sea'] = (AS_condition).sum().values
region_cell_counts['Bay of Bengal'] = (BB_condition).sum().values
region_cell_counts['Indian Ocean'] = (INDO_condition).sum().values
region_cell_counts['Southeastern Asia'] = (SEAS_condition).sum().values
region_cell_counts['North Pacific'] = (NPAC_condition).sum().values
region_cell_counts['Arctic'] = (ARCT_condition).sum().values
region_cell_counts['Australia/South Pacific'] = (AUSP_condition).sum().values
region_cell_counts['Southern Ocean'] = (SO_condition).sum().values
region_cell_counts['Central Pacific / Asia'] = (CPAO_condition).sum().values

total_cell_counts['total'] = all_FE_DEP.dims['lat'] * all_FE_DEP.dims['lon']

# Calculate the total number of grid cells across all regions
total_grid_cells_regional_sums = sum(region_cell_counts.values()) - 5248
total_grid_cells = total_cell_counts['total'] 

# --------- EMISSIONS BUDGETS
# Loop over all variables in the dataset
for var_name in all_FE_DEP.data_vars:
    if "DEP" in var_name:  
        # Calculate individual emissions (for each cell) and convert from sec to annual and from kg to Tg
        individual_cell_emissions[var_name] = (all_FE_DEP[var_name] * all_FE_DEP['cell_area'] * 3600 * 24 * 365 *1E-9)
        
       # Calculate the total budget by summing individual emissions
        total_budget = float(individual_cell_emissions[var_name].sum().values)  # Ensure it's a float

        # Calculate Regional budgets by summing specified grid cells
        SATL_budget = float(individual_cell_emissions[var_name].where(SATL_condition).sum().values)
        NATL_budget = float(individual_cell_emissions[var_name].where(NATL_condition).sum().values)
        AS_budget = float(individual_cell_emissions[var_name].where(AS_condition).sum().values)
        BB_budget = float(individual_cell_emissions[var_name].where(BB_condition).sum().values)
        INDO_budget = float(individual_cell_emissions[var_name].where(INDO_condition).sum().values)
        SEAS_budget = float(individual_cell_emissions[var_name].where(SEAS_condition).sum().values)
        NPAC_budget = float(individual_cell_emissions[var_name].where(NPAC_condition).sum().values)
        ENPAC_budget = float(individual_cell_emissions[var_name].where(NPAC_condition).sum().values)
        WNPAC_budget = float(individual_cell_emissions[var_name].where(NPAC_condition).sum().values)
        ARCT_budget = float(individual_cell_emissions[var_name].where(ARCT_condition).sum().values)
        AUSP_budget = float(individual_cell_emissions[var_name].where(AUSP_condition).sum().values)
      #  SIND_budget = float(individual_cell_emissions[var_name].where(SIND_condition).sum().values)
        SO_budget = float(individual_cell_emissions[var_name].where(SO_condition).sum().values)
        CPAO_budget = float(individual_cell_emissions[var_name].where(CPAO_condition).sum().values)
        
        # Store the budgets for the variable in a nested dictionary
        global_budgets[var_name] = {
            "Total_Budget": total_budget,
            "SATL_Budget": SATL_budget,
            "NATL_Budget": NATL_budget,
            "AS_Budget": AS_budget,
            "BB_Budget": BB_budget,
            "INDO_budget": INDO_budget,
            "SEAS_Budget": SEAS_budget,
            "NPAC_Budget": NPAC_budget,
            "ENPAC_Budget": ENPAC_budget,
            "WNPAC_Budget": WNPAC_budget,
            "ARCT_Budget": ARCT_budget,
           "AUSP_Budget": AUSP_budget,
      #     "SIND_Budget": SIND_budget,
            "SO_Budget": SO_budget,
            "CPAO_Budget": CPAO_budget
        }

budget_df = pd.DataFrame(global_budgets).T  # Transpose to get variables as rows and budgets as columns

# Reset index and give the proper column name for variables
budget_df.reset_index(inplace=True)
budget_df.rename(columns={'index': 'Variable'}, inplace=True)

# Ensure all columns have the correct numeric type
budget_df = budget_df.apply(pd.to_numeric, errors='ignore')  # Apply to convert all numeri

# Save the total budgets and individual emissions to separate sheets in an Excel file
output_file = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\dep_budgets_PD_V4.xlsx"

with pd.ExcelWriter(output_file) as writer:
   budget_df.to_excel(writer, sheet_name="Depo_Budgets", index=False)  
   print(f"Deposition budgets saved to {output_file}")

#%% --- AND NOW FOR EMISSIONS -------------------------------------------------------------------------------
#ds_A= xr.open_dataset("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Transient Run\\input data files\\ClaquinMineralsCAM6_SPEWFUELS-Dec2023.nc")

#ds_A= xr.open_dataset("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Transient Run\\input data files\\SPEWFUEL_emissions_1850-2010_current_CF.nc") 

#ds_A= xr.open_dataset("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Transient Run\\input data files\\SPEWFUEL_emissions_1850-2010_CF_no_BF_globally.nc")

ds_A= xr.open_dataset("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Transient Run\\input data files\\SPEWFUEL_emissions_1850-2010_CF_no_BF_in_SouthernHemi_5.nc")

selected_vars = ["Fine_InsolFe_Comb", "Coar_InsolFe_Comb", "Fine_SolFe_Comb", "Coar_SolFe_Comb"]  
all_FE_DEP = ds_A[selected_vars]

all_FE_DEP = all_FE_DEP.sel(time="2010")
all_FE_DEP = all_FE_DEP.mean(dim="time")

# FIND AREA OF INDIVIDUAL GRID CELLS --------------------------------
lon = all_FE_DEP['lon'].values  # Longitude in degrees
lat = all_FE_DEP['lat'].values  # Latitude in degrees

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

# add cell area to original xarrays for easy calling 
all_FE_DEP['cell_area'] = xr.DataArray(
    cell_areas_staggered,
    dims=["lat", "lon"],  # Same dimensions as in the original dataset
    coords={"lat": all_FE_DEP['lat'], "lon": all_FE_DEP['lon']},  # Use original coordinates
    attrs={
        "units": "m^2",  # Specify units for the cell area
        "description": "Calculated grid cell area using staggered grid approach",
    },
)

# ---- ASSIGNING ATM. DEPOSITION REGIONS 
import pandas as pd
# Dictionaries to store both individual and total budgets
global_budgets = {}
individual_cell_emissions = {}

# Define regional conditions
# DO NOT CHANGE THESE!!! THEY ARE GOOD NOW FOR ADDING FINAL BUDGETS TO TOTAL 
# South Atlantic
#SATL_condition = (all_FE_DEP['lon'] >= 295.0) & (all_FE_DEP['lat'] <= -15.0) & (all_FE_DEP['lat'] > -48.0) | \
#                 (all_FE_DEP['lon'] < 30) & (all_FE_DEP['lat'] <= -15.0) & (all_FE_DEP['lat'] >-48.0)
# North Atlantic
#NATL_condition = (all_FE_DEP['lon'] > 265.0) & (all_FE_DEP['lat'] > 29.0) & (all_FE_DEP['lat'] < 60.0) | \
#                 (all_FE_DEP['lon'] < 110.0) & (all_FE_DEP['lat'] > 29.0) & (all_FE_DEP['lat'] < 60.0)
# Arabian Sea
#AS_condition = (all_FE_DEP['lon'] < 78.0) & (all_FE_DEP['lon'] >= 30.0) & \
 #              (all_FE_DEP['lat'] <= 29.0) & (all_FE_DEP['lat'] > -48.0)
# Bay of Bengal
#BB_condition = (all_FE_DEP['lon'] < 110.0) & (all_FE_DEP['lon'] >= 78.0) & \
#               (all_FE_DEP['lat'] <= 29.0) & (all_FE_DEP['lat'] > -48.0)
# Indian Ocean 
#INDO_condition = (all_FE_DEP['lon'] < 110.0) & (all_FE_DEP['lon'] >= 30.0) & \
#               (all_FE_DEP['lat'] <= 29.0) & (all_FE_DEP['lat'] > -48.0)
# Southeastern Asia
#SEAS_condition = (all_FE_DEP['lon'] < 150.0) & (all_FE_DEP['lon'] >= 110.0) & \
 #              (all_FE_DEP['lat'] < 60.0) & (all_FE_DEP['lat'] > -15.0)
## North Pacific
#NPAC_condition = (all_FE_DEP['lon'] <= 265.0) & (all_FE_DEP['lon'] >= 150.0) & \
#               (all_FE_DEP['lat'] < 60.0) & (all_FE_DEP['lat'] > 29.0)
# North Pacific
#ENPAC_condition = (all_FE_DEP['lon'] <= 180.0) & (all_FE_DEP['lon'] >= 150.0) & \
 #              (all_FE_DEP['lat'] < 60.0) & (all_FE_DEP['lat'] > 29.0)
# North Pacific
#WNPAC_condition = (all_FE_DEP['lon'] <= 265.0) & (all_FE_DEP['lon'] >= 180.0) & \
#               (all_FE_DEP['lat'] < 60.0) & (all_FE_DEP['lat'] > 29.0)
# Arctic
#ARCT_condition = (all_FE_DEP['lon'] >= 0.0) & (all_FE_DEP['lat'] >= 60.0) &\
#               (all_FE_DEP['lon'] >= 0.0)  & (all_FE_DEP['lat'] <= 90.0)
# Australia/South Pacific
#AUSP_condition = (all_FE_DEP['lon'] < 295.0) & (all_FE_DEP['lon'] >= 110.0) & \
#               (all_FE_DEP['lat'] <= -15.0) & (all_FE_DEP['lat'] > -48.0)
# Southern Ocean
#SO_condition = (all_FE_DEP['lon'] >= 0.0) & (all_FE_DEP['lat'] <= -48.0) &\
 #              (all_FE_DEP['lon'] >= 0.0)  & (all_FE_DEP['lat'] >= -90.0)
# Central Pacific / Asia
#CPAO_condition = (all_FE_DEP['lon'] >= 150.0) & (all_FE_DEP['lat'] > -15.0) & (all_FE_DEP['lat'] <= 29.0) | \
#                 (all_FE_DEP['lon'] < 30.0) & (all_FE_DEP['lat'] > -15.0) & (all_FE_DEP['lat'] <= 29.0)

# Continental emission settings
# Defining Regional Boundaries (mostly separated by continents)
# South America
SAM_condition = (all_FE_DEP['lon'] >= 279.0) & (all_FE_DEP['lon'] <= 326.0) & \
                (all_FE_DEP['lat'] >= -56.0) & (all_FE_DEP['lat'] <= 12.0)

# North America other than USA
NAM_condition = ( # MEXICO / CENTRAL AM
                ((all_FE_DEP['lon'] >= 190.0) & (all_FE_DEP['lon'] <= 310.0) & \
                (all_FE_DEP['lat'] > 12.0) & (all_FE_DEP['lat'] < 30.0))
                | # CANADA
                ((all_FE_DEP['lon'] >= 190.0) & (all_FE_DEP['lon'] <= 310.0) & \
                (all_FE_DEP['lat'] > 49.0) & (all_FE_DEP['lat'] < 60.0))
                )       

USA_condition = (all_FE_DEP['lon'] >= 190.0) & (all_FE_DEP['lon'] <= 310.0) & \
                (all_FE_DEP['lat'] >= 30.0) & (all_FE_DEP['lat'] <= 49.0)

# Africa -- when straddling prime meridian must code this way 
AFR_condition = (
                ((all_FE_DEP['lon'] >= 340) & (all_FE_DEP['lat'] < 37.0) & (all_FE_DEP['lat'] >= 0.0) | \
                (all_FE_DEP['lon'] <= 35) & (all_FE_DEP['lat'] < 37.0) & (all_FE_DEP['lat'] >= 0.0))
                |
                ((all_FE_DEP['lon'] >= 340) & (all_FE_DEP['lat'] < 15.0) & (all_FE_DEP['lat'] >= 0.0) | \
                (all_FE_DEP['lon'] <= 50) & (all_FE_DEP['lat'] < 15.0) & (all_FE_DEP['lat'] >= 0.0))
                )

# Southern portions of Africa < Equator
SAFR_condition = (all_FE_DEP['lon'] >= 340) & (all_FE_DEP['lat'] < 0.0) & (all_FE_DEP['lat'] >= -35.0) | \
                (all_FE_DEP['lon'] <= 50) & (all_FE_DEP['lat'] < 0.0) & (all_FE_DEP['lat'] >= -35.0)

# Europe 
EUR_condition = (all_FE_DEP['lon'] >= 335) & (all_FE_DEP['lat'] < 60.0) & (all_FE_DEP['lat'] >= 37.0) | \
                (all_FE_DEP['lon'] <= 35) & (all_FE_DEP['lat'] < 60.0) & (all_FE_DEP['lat'] >= 37.0)

# Northern/ Eastern Asia/China, Himalayas separate the air mass movement so I need to separate finer around them
CHINA_condition = (
                  ((all_FE_DEP['lon'] < 180) & (all_FE_DEP['lon'] > 70.0) & \
                  (all_FE_DEP['lat'] < 55.0) & (all_FE_DEP['lat'] >= 30))
                  | 
                  ((all_FE_DEP['lon'] < 180) & (all_FE_DEP['lon'] > 85.0) & \
                  (all_FE_DEP['lat'] < 30.0) & (all_FE_DEP['lat'] >= 22.0))
                  )   

# Southeastern Asia and India
SEAS_condition =  (
                  ((all_FE_DEP['lon'] < 180) & (all_FE_DEP['lon'] > 55.0) & \
                  (all_FE_DEP['lat'] < 30) & (all_FE_DEP['lat'] >= 22.0))
                  | 
                  ((all_FE_DEP['lon'] < 180) & (all_FE_DEP['lon'] > 55.0) & \
                  (all_FE_DEP['lat'] < 22.0) & (all_FE_DEP['lat'] >= -10.0))
                  )

# Australia/South Pacific
AUS_condition = (all_FE_DEP['lon'] <= 180) & (all_FE_DEP['lon'] >= 110.0) & \
               (all_FE_DEP['lat'] < -10.0) & (all_FE_DEP['lat'] > -50.0)

# Southern Hemisphere
SH_condition = (all_FE_DEP['lon'] <= 360) & \
               (all_FE_DEP['lat'] < 0.0) & (all_FE_DEP['lat'] >= -90.0)

# Initialize an empty dictionary to store the grid cell counts
region_cell_counts = {}
total_cell_counts = {}

# Loop through each region condition and count the grid cells that satisfy the condition
#region_cell_counts['South Atlantic'] = (SATL_condition).sum().values
#region_cell_counts['North Atlantic'] = (NATL_condition).sum().values
#region_cell_counts['Arabian Sea'] = (AS_condition).sum().values
#region_cell_counts['Bay of Bengal'] = (BB_condition).sum().values
#region_cell_counts['Indian Ocean'] = (INDO_condition).sum().values
#region_cell_counts['Southeastern Asia'] = (SEAS_condition).sum().values
#region_cell_counts['North Pacific'] = (NPAC_condition).sum().values
#region_cell_counts['ENorth Pacific'] = (ENPAC_condition).sum().values
#region_cell_counts['WNorth Pacific'] = (WNPAC_condition).sum().values
#region_cell_counts['Arctic'] = (ARCT_condition).sum().values
#region_cell_counts['Australia/South Pacific'] = (AUSP_condition).sum().values
#region_cell_counts['Southern Ocean'] = (SO_condition).sum().values
#region_cell_counts['Central Pacific / Asia'] = (CPAO_condition).sum().values

region_cell_counts['South America'] = (SAM_condition).sum().values
region_cell_counts['North America'] = (NAM_condition).sum().values
region_cell_counts['USA'] = (USA_condition).sum().values
region_cell_counts['Africa'] = (AFR_condition).sum().values
region_cell_counts['Southern Africa'] = (SAFR_condition).sum().values
region_cell_counts['Europe'] = (EUR_condition).sum().values
region_cell_counts['China'] = (CHINA_condition).sum().values
region_cell_counts['SE Asia'] = (SEAS_condition).sum().values
region_cell_counts['Australia'] = (AUS_condition).sum().values

total_cell_counts['total'] = all_FE_DEP.dims['lat'] * all_FE_DEP.dims['lon']

# Calculate the total number of grid cells across all regions
total_grid_cells_regional_sums = sum(region_cell_counts.values()) - 5248
total_grid_cells = total_cell_counts['total'] 

# --------- EMISSIONS BUDGETS
# Loop over all variables in the dataset
for var_name in all_FE_DEP.data_vars:
    if "Fe" in var_name:  
        # Calculate individual emissions (for each cell) and convert from sec to annual and from kg to Tg
        individual_cell_emissions[var_name] = (all_FE_DEP[var_name] * all_FE_DEP['cell_area'] * 3600 * 24 * 365 *1E-9)
        
       # Calculate the total budget by summing individual emissions
        total_budget = float(individual_cell_emissions[var_name].sum().values)  # Ensure it's a float

        # Calculate Regional budgets by summing specified grid cells
        SAM_budget = float(individual_cell_emissions[var_name].where(SAM_condition).sum().values)
        NAM_budget = float(individual_cell_emissions[var_name].where(NAM_condition).sum().values)
        USA_budget = float(individual_cell_emissions[var_name].where(USA_condition).sum().values)
        AFR_budget = float(individual_cell_emissions[var_name].where(AFR_condition).sum().values)
        SAFR_budget = float(individual_cell_emissions[var_name].where(SAFR_condition).sum().values)
        EUR_budget = float(individual_cell_emissions[var_name].where(EUR_condition).sum().values)
        CHINA_budget = float(individual_cell_emissions[var_name].where(CHINA_condition).sum().values)
        SEAS_budget = float(individual_cell_emissions[var_name].where(SEAS_condition).sum().values)
        AUS_budget = float(individual_cell_emissions[var_name].where(AUS_condition).sum().values)
        
        # Store the budgets for the variable in a nested dictionary
        global_budgets[var_name] = {
            "Total_Budget": total_budget,
            "SAM_budget": SAM_budget,
            "NAM_budget": NAM_budget,
            "USA_budget": USA_budget,
            "AFR_budget": AFR_budget,
            "SAFR_budget": SAFR_budget,
            "EUR_budget": EUR_budget,
            "CHINA_budget": CHINA_budget,
            "SEAS_budget": SEAS_budget,
            "AUS_budget": AUS_budget,
        }

budget_df = pd.DataFrame(global_budgets).T  # Transpose to get variables as rows and budgets as columns

# Reset index and give the proper column name for variables
budget_df.reset_index(inplace=True)
budget_df.rename(columns={'index': 'Variable'}, inplace=True)

# Ensure all columns have the correct numeric type
budget_df = budget_df.apply(pd.to_numeric, errors='ignore')  # Apply to convert all numeri

# Save the total budgets and individual emissions to separate sheets in an Excel file
output_file = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Transient Run\\input data files\\EMISS_budgets_NA_CF_test.xlsx"

with pd.ExcelWriter(output_file) as writer:
   budget_df.to_excel(writer, sheet_name="Emiss_Budgets", index=False)  
   print(f"Emission budgets saved to {output_file}")
#%% SIMPLIFIED EMISSION BUDGETS -- For regional and global budgets 
import numpy as np
import xarray as xr
from netCDF4 import Dataset

# LOAD IN THE NETCDF FILE -----------------------------------------
#sim_directory = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\ClaquinMineralsCAM6-LL_SPEW_BGC_Feb2020-SDR_PI_bf.nc"
#ds_h1 = xr.open_dataset(sim_directory)

# EXTRACT VARIABLES ---------------------------------------------------------------
# IRON PI 
#FE_EMISS_COA = ds_h1["FFBFc"]
#FE_EMISS_FINE = ds_h1["FFBFc"]

# IRON PD AND FU 
#CoaMed_FCOALFF = ds_h1["CoaMed_FCOALFF"]
#FineMed_FCOALFF = ds_h1["CoaMed_FCOALFF"]

#CoaMed_FCOALBF = ds_h1["CoaMed_FCOALBF"]
#FineMed_FCOALBF = ds_h1["CoaMed_FCOALBF"]

#CoaMed_FWOOD = ds_h1["CoaMed_FWOOD"]
#FineMed_FWOOD = ds_h1["CoaMed_FWOOD"]

#CoaMed_FOIL = ds_h1["CoaMed_FOIL"]
#FineMed_FOIL = ds_h1["CoaMed_FOIL"]

#CoaMed_FSMELT = ds_h1["CoaMed_FSMELT"]
#FineMed_FSMELT = ds_h1["CoaMed_FSMELT"]

#FE_EMISS_COA = CoaMed_FCOALFF + CoaMed_FCOALBF + CoaMed_FWOOD + CoaMed_FOIL + CoaMed_FSMELT
#FE_EMISS_FINE = FineMed_FCOALFF + FineMed_FCOALBF + FineMed_FWOOD + FineMed_FOIL + FineMed_FSMELT

# BC is more complicated 
import xarray as xr
import pandas as pd
import numpy as np

# Load the NetCDF file -- BC EMISSIONS SSP370
# ds_SSP370 = xr.open_dataset("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\BC-em-anthro_input4MIPs_emissions_ScenarioMIP_IAMC-AIM-ssp370-1-1_gn_2090_2100MEAN_remapcon_regridded.nc")
ds_SSP370 = xr.open_dataset("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\BC-em-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_2000-2015MEAN_remapcon_regridded.nc") 
ds_SSP370 = ds_SSP370.drop_vars(["time_bnds", "sector_bnds", "time"])
ds_SSP370 = ds_SSP370.squeeze(dim='time') 
# # BOUNDING HAS BEEN AN ISSUE -- despite remapping, final value was 90.00000058 instead of 90.

# Load the NetCDF file -- BC EMISSIONS PD CMIP6
#ds_PD = xr.open_dataset("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\BC-em-anthro_input4MIPs_emissions_CMIP_CEDS-2017-05-18_gn_2000-2015MEAN_remapcon_regridded.nc") 
#ds_PD = ds_PD.drop_vars(["time_bnds", "sector_bnds", "time"])
#ds_PD = ds_PD.drop_dims(["time"])
#ds_PD = ds_PD.squeeze(dim='time') 

# Load the NetCDF file -- Fe EMISSIONS PD CMIP6 (use as OCNFRAC too)
ds_PD_Fe = xr.open_dataset("C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\ClaquinMineralsCAM6_SPEWFUELS-Dec2023_OCNFRAC_remapcon_regridded.nc")
ds_PD_Fe['OCNFRAC'] = ds_PD_Fe['OCNFRAC'].squeeze(dim='time')
ds_PD_Fe = ds_PD_Fe.drop_dims(["time"]) 

# Extract Black Carbon (BC) emission values SSP370 
BC_em_anthro_SSP370 = ds_SSP370['BC_em_anthro'] # Emission data
lon = ds_SSP370['lon'].values  # Longitude values -- the .values extracts as a numPy array, removes metadata
lat = ds_SSP370['lat'].values  # Latitude values
sector = ds_SSP370['sector'].values # Sector identifiers

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

BC_COALBF_SSP370 = BC_COALBF_SSP370.fillna(0.0)
BC_WOOD_SSP370 = BC_WOOD_SSP370.fillna(0.0)
BC_COALFF_SSP370 = (BC_emiss_1_SSP370 + BC_emiss_2_SSP370)
BC_OIL_SSP370 = filtered_BC_emiss_3_SSP370.fillna(filtered_BC_emiss_7_SSP370)

BC = BC_COALBF_SSP370 + BC_WOOD_SSP370 + BC_COALFF_SSP370 + BC_OIL_SSP370

# ADD ALL VARIABLES TO SINGLE DATASET --------------------------------------------
all_FE_DEP = xr.Dataset({
    # TOTAL SOURCES
   # "FE_EMISS_COA": FE_EMISS_COA,
   # "FE_EMISS_FINE": FE_EMISS_FINE,
    "BC": BC,
    })

# FIND AREA OF INDIVIDUAL GRID CELLS --------------------------------
lon = all_FE_DEP['lon'].values  # Longitude in degrees
lat = all_FE_DEP['lat'].values  # Latitude in degrees

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

# add cell area to original xarrays for easy calling 
all_FE_DEP['cell_area'] = xr.DataArray(
    cell_areas_staggered,
    dims=["lat", "lon"],  # Same dimensions as in the original dataset
    coords={"lat": all_FE_DEP['lat'], "lon": all_FE_DEP['lon']},  # Use original coordinates
    attrs={
        "units": "m^2",  # Specify units for the cell area
        "description": "Calculated grid cell area using staggered grid approach",
    },
)

# ---- ASSIGNING ATM. DEPOSITION REGIONS 
import pandas as pd
# Dictionaries to store both individual and total budgets
global_budgets = {}
individual_cell_emissions = {}

# Define regional conditions
# DO NOT CHANGE THESE!!! THEY ARE GOOD NOW FOR ADDING FINAL BUDGETS TO TOTAL 
# South Atlantic
SATL_condition = (all_FE_DEP['lon'] >= 295.0) & (all_FE_DEP['lat'] <= -15.0) & (all_FE_DEP['lat'] > -48.0) | \
                 (all_FE_DEP['lon'] < 30) & (all_FE_DEP['lat'] <= -15.0) & (all_FE_DEP['lat'] >-48.0)
# North Atlantic
NATL_condition = (all_FE_DEP['lon'] > 265.0) & (all_FE_DEP['lat'] > 29.0) & (all_FE_DEP['lat'] < 60.0) | \
                 (all_FE_DEP['lon'] < 110.0) & (all_FE_DEP['lat'] > 29.0) & (all_FE_DEP['lat'] < 60.0)
# Arabian Sea
AS_condition = (all_FE_DEP['lon'] < 78.0) & (all_FE_DEP['lon'] >= 30.0) & \
               (all_FE_DEP['lat'] <= 29.0) & (all_FE_DEP['lat'] > -48.0)
# Bay of Bengal
BB_condition = (all_FE_DEP['lon'] < 110.0) & (all_FE_DEP['lon'] >= 78.0) & \
               (all_FE_DEP['lat'] <= 29.0) & (all_FE_DEP['lat'] > -48.0)
# Indian Ocean 
INDO_condition = (all_FE_DEP['lon'] < 110.0) & (all_FE_DEP['lon'] >= 30.0) & \
               (all_FE_DEP['lat'] <= 29.0) & (all_FE_DEP['lat'] > -48.0)
# Southeastern Asia
SEAS_condition = (all_FE_DEP['lon'] < 150.0) & (all_FE_DEP['lon'] >= 110.0) & \
               (all_FE_DEP['lat'] < 60.0) & (all_FE_DEP['lat'] > -15.0)
# North Pacific
NPAC_condition = (all_FE_DEP['lon'] <= 265.0) & (all_FE_DEP['lon'] >= 150.0) & \
               (all_FE_DEP['lat'] < 60.0) & (all_FE_DEP['lat'] > 29.0)
# Arctic
ARCT_condition = (all_FE_DEP['lon'] >= 0.0) & (all_FE_DEP['lat'] >= 60.0) &\
               (all_FE_DEP['lon'] >= 0.0)  & (all_FE_DEP['lat'] <= 90.0)
# Australia/South Pacific
AUSP_condition = (all_FE_DEP['lon'] < 295.0) & (all_FE_DEP['lon'] >= 110.0) & \
               (all_FE_DEP['lat'] <= -15.0) & (all_FE_DEP['lat'] > -48.0)
# Southern Ocean
SO_condition = (all_FE_DEP['lon'] >= 0.0) & (all_FE_DEP['lat'] <= -48.0) &\
               (all_FE_DEP['lon'] >= 0.0)  & (all_FE_DEP['lat'] >= -90.0)
# Central Pacific / Asia
CPAO_condition = (all_FE_DEP['lon'] >= 150.0) & (all_FE_DEP['lat'] > -15.0) & (all_FE_DEP['lat'] <= 29.0) | \
                 (all_FE_DEP['lon'] < 30.0) & (all_FE_DEP['lat'] > -15.0) & (all_FE_DEP['lat'] <= 29.0)

# Initialize an empty dictionary to store the grid cell counts
region_cell_counts = {}
total_cell_counts = {}

# Loop through each region condition and count the grid cells that satisfy the condition
region_cell_counts['South Atlantic'] = (SATL_condition).sum().values
region_cell_counts['North Atlantic'] = (NATL_condition).sum().values
region_cell_counts['Arabian Sea'] = (AS_condition).sum().values
region_cell_counts['Bay of Bengal'] = (BB_condition).sum().values
region_cell_counts['Indian Ocean'] = (INDO_condition).sum().values
region_cell_counts['Southeastern Asia'] = (SEAS_condition).sum().values
region_cell_counts['North Pacific'] = (NPAC_condition).sum().values
region_cell_counts['Arctic'] = (ARCT_condition).sum().values
region_cell_counts['Australia/South Pacific'] = (AUSP_condition).sum().values
region_cell_counts['Southern Ocean'] = (SO_condition).sum().values
region_cell_counts['Central Pacific / Asia'] = (CPAO_condition).sum().values

total_cell_counts['total'] = all_FE_DEP.dims['lat'] * all_FE_DEP.dims['lon']

# Calculate the total number of grid cells across all regions
total_grid_cells_regional_sums = sum(region_cell_counts.values()) - 5248
total_grid_cells = total_cell_counts['total'] 

# --------- EMISSIONS BUDGETS
# Loop over all variables in the dataset
for var_name in all_FE_DEP.data_vars:
    if "BC" in var_name:  
        # Calculate individual emissions (for each cell) and convert from sec to annual and from kg to Tg
        individual_cell_emissions[var_name] = (all_FE_DEP[var_name] * all_FE_DEP['cell_area'] * 3600 * 24 * 365 *1E-9)
        
       # Calculate the total budget by summing individual emissions
        total_budget = float(individual_cell_emissions[var_name].sum().values)  # Ensure it's a float

        # Calculate Regional budgets by summing specified grid cells
        SATL_budget = float(individual_cell_emissions[var_name].where(SATL_condition).sum().values)
        NATL_budget = float(individual_cell_emissions[var_name].where(NATL_condition).sum().values)
        AS_budget = float(individual_cell_emissions[var_name].where(AS_condition).sum().values)
        BB_budget = float(individual_cell_emissions[var_name].where(BB_condition).sum().values)
        INDO_budget = float(individual_cell_emissions[var_name].where(INDO_condition).sum().values)
        SEAS_budget = float(individual_cell_emissions[var_name].where(SEAS_condition).sum().values)
        NPAC_budget = float(individual_cell_emissions[var_name].where(NPAC_condition).sum().values)
        ARCT_budget = float(individual_cell_emissions[var_name].where(ARCT_condition).sum().values)
        AUSP_budget = float(individual_cell_emissions[var_name].where(AUSP_condition).sum().values)
       # SIND_budget = float(individual_cell_emissions[var_name].where(SIND_condition).sum().values)
        SO_budget = float(individual_cell_emissions[var_name].where(SO_condition).sum().values)
        CPAO_budget = float(individual_cell_emissions[var_name].where(CPAO_condition).sum().values)
        
        # Store the budgets for the variable in a nested dictionary
        global_budgets[var_name] = {
            "Total_Budget": total_budget,
            "SATL_Budget": SATL_budget,
            "NATL_Budget": NATL_budget,
            "AS_Budget": AS_budget,
            "BB_Budget": BB_budget,
            "INDO_budget": INDO_budget,
            "SEAS_Budget": SEAS_budget,
            "NPAC_Budget": NPAC_budget,
            "ARCT_Budget": ARCT_budget,
            "AUSP_Budget": AUSP_budget,
          #  "SIND_Budget": SIND_budget,
            "SO_Budget": SO_budget,
            "CPAO_Budget": CPAO_budget
        }

budget_df = pd.DataFrame(global_budgets).T  # Transpose to get variables as rows and budgets as columns

# Reset index and give the proper column name for variables
budget_df.reset_index(inplace=True)
budget_df.rename(columns={'index': 'Variable'}, inplace=True)

# Ensure all columns have the correct numeric type
budget_df = budget_df.apply(pd.to_numeric, errors='ignore')  # Apply to convert all numeri

# Save the total budgets and individual emissions to separate sheets in an Excel file
output_file = "C:\\Users\\heplaas\\OneDrive - North Carolina State University\\Collaborations\\Coal Fly Ash\\data\\PD_BC_EMISS_budgets.xlsx"

with pd.ExcelWriter(output_file) as writer:
   budget_df.to_excel(writer, sheet_name="Emiss_Budgets", index=False)  
   print(f"Emission budgets saved to {output_file}")

