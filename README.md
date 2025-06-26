# CombustionFe

R and Python scripts as utilized for the publication entitled: 
Residential coal and biofuel burning is a signficant source of 
soluble aerosol iron to the atmosphere, in preparation for 
Atmospheric Chemistry and Physics - EGU. Supporting data files 
are provided in this repository.

All data and programming scripts are labeled with 
Creative Commons Attribution 4.0 labeled for cited reuse.

Data files include Fe and BC emissions inventories for the 
pre-industrial (PI), present day (PD), and future (FU) eras, 
as well as total and soluble Fe deposition fluxes for the 
same emissions scenarios. Model input/output data is provided 
in a 3D lat-lon grid encoded in a netcdf file. Cleaned model 
deposition fluxes and surface concentrations for Fe, with 
observations of soluble Fe used in model-observation comparison 
are provided in excel file format. 

Model output file naming conventions (made available upon request
due to storage space requirements -- please contact Haley E. Plaas
at hep2126@columbia.edu and cc haleyeplaas@gmail.com): 
CAM6-MIMI-a-INDCOALb-RESICOALc-WOODd-OILe.cam.hf.2009-2011.nc ; 
where 

a = emissions scenario PI , PD, or FU (SSP370)

b / c / d / e = fractional Fe solubility applied and represented 
as a percentage for the specified fuel type, with 
INDCOAL = industrial coal, 
RESICOAL = residential coal, 
WOOD = residential biofuel burning, and 
OIL = oil (shipping + transportation)

f = 1 / 2 , with 1 representing source apportioned dry and wet 
deposition fluxes of total and soluble Fe and 2 representing 
source apportioned total and soluble Fe surface concentrations
