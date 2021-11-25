# noise_raster

This repository contains python code to:

1. create a list of rasters contained within each source folder path
2. validate source raster data: file extension check, coordinate reference system check, raster extent check, convert negative and 0 values to no data and logs negative pixels lower than -90, set all no data values to -99, set all data types to 32 bit
3. reproject source raster data to EPSG:3035
4. merge each list of source rasters to a single merged raster for each noise type: resample source raster data with different resolutions to the highest resolution
5. create a multi-band virtual raster of all noise sources
6. create a masked array based on the multi-band raster to handle no data values
7. energetically add raster data from different noise sources
8. write added array to new raster
9. reclassify and vectorize added raster data


