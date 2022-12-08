# USGS Streamflow Data

This folder contains code to process and download streamgage data
available from the National Water Information Service (NWIS) at the
United States Geological Service (USGS). Flow data was retrieved from
the NWIS [1] using the R package *r**n**o**a**a* [2]. Drainage area
geometries were retrieved from Streamstats [3] using both the API and
the batch processing tool. The file `get_streamflow.Rmd` guides the user
through the process of downloading flows, downloading drainage area
geometries, converting streamflow to runoff, and regionalizing runoff to
match the resolution of the MERRA-2/SPEAR 50km × 50km grid using the
method from Brakebill et al. (2011) [4].

<br><br>

[1] USGS. (2022). National Water Information System. Retrieved from
<https://waterdata.usgs.gov/nwis>.

[2] Chamberlain, S. (2021). rnoaa: “NOAA” Weather Data from R. Retrieved
from <https://cran.r-project.org/package=rnoaa>.

[3] USGS. (2019). StreamStats. Retrieved from
<http://streamstats.usgs.gov/>.

[4] Brakebill, J. W., Wolock, D. M., & Terziotti, S. E. (2011). Digital
Hydrologic Networks Supporting Applications Related to Spatially
Referenced Regression Modeling. Journal of the American Water Resources
Association, 47(5), 916–932.
<https://doi.org/10.1111/j.1752-1688.2011.00578.x>.
