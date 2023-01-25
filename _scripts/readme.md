# Scripts Folder

This folder contains three data processing scripts.

-   `create_df_merra.Rmd` creates lists of dataframes at both a 3-hour
    (**df\_3hr**) and a 24-hour (**df\_24hr**) resolution for all
    MERRA-2 grid cells within California, where each grid cell is
    indexed to one dataframe containing timeseries values from
    1981-2022.
-   `create_df_gfdl.Rmd` creates three lists of dataframes. **df\_hist**
    is a list of dataframes for all GFDL SPEAR grid cells within
    California, where each grid cell is indexed to one dataframe
    containing timeseries values from 1921-2010, and each dataframe in
    the list contains daily data from five GFDL SPEAR ensemble members
    under historical climate forcings. **df\_ssp245** and **df\_ssp585**
    are lists of datasets constructed in the same fashion for SSP 2–4.5
    and 5–8.5 forcings, respectively, and cover the time period from
    2021-2090.
-   `create_df_functions.R` is a helper script containing functions
    common to both `.Rmd` scripts.

All of the datasets created in `create_df_gfdl.Rmd` use bias-corrected
IVT data as presented in the `_data/GFDL/biascorr` folder. GFDL SPEAR
data is recorded at a 24-hour resolution.

All datasets contain the following variables: IVT, 5-day rolling IVT, AR
indicator, AR ID number, AR maximum IVT, AR duration, sequence
indicator, sequence ID number, sequence maximum IVT, and sequence
duration. In addition, the datasets created in `create_df_merra.Rmd`
contain variables for precipitation (mm) and soil moisture (mm/m).
