# noise_raster

This repository contains python code to:

1. create a list of rasters contained within each source folder path
2. validate source raster data: file extension check, coordinate reference system check, raster extent check, set all no data values to -99, set all data types to 32 bit
3. merge each list of source rasters to a single merged raster for each noise type: resample source raster data with different resolutions to the highest resolution
4. create a multi-band virtual raster of all noise sources when there are more than one noise source
5. create an array based on the multi-band raster and convert all negative values to 0 as input to addition
6. energetically add raster data from different noise sources
7. write added array to new raster
8. reclassify and vectorize added raster data
9. reproject raster to EPSG:3035


