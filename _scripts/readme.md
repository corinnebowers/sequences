# Scripts Folder

This folder contains three data processing scripts.

`create_df_merra.Rmd` creates a list of dataframes at both a 3-hour
(*df\_3hr*) and a 24-hour (*df\_24hr*) resolution for all MERRA-2 grid
cells within California, where each grid cell is indexed to one
dataframe containing timeseries values from 1981-2022.

`create_df_gfdl.Rmd` creates three lists of dataframes. *df\_hist*
contains daily data from five GFDL SPEAR ensemble members under
historical climate forcings. It is a list of dataframes for all GFDL
SPEAR grid cells within California, where each grid cell is indexed to
one dataframe containing timeseries values from 1921-2010.  
*df\_ssp245* and *df\_ssp585* are the same for SSP 2–4.5 and 5–8.5
forcings, respectively, and cover the time period from 2021-2090.

All datasets contain the following variables:

-   IVT (kg/m/s)
-   5-day rolling IVT
-   AR indicator
-   AR ID number
-   AR maximum IVT
-   AR duration
-   Sequence indicator
-   Sequence ID number
-   Sequence maximum IVT
-   Sequence duration

All of the datasets created in `create_df_gfdl.Rmd` use bias-corrected
IVT data as presented in the `_data/GFDL/biascorr` folder. GFDL SPEAR
data is recorded at a 24-hour resolution. In addition to the variables
above, the datasets created in `create_df_merra.Rmd` contain variables
for precipitation (mm) and soil moisture (mm/m).
