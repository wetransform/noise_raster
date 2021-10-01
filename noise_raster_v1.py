# -*- coding: utf-8 -*-
"""
/***************************************************************************
 NoiseRaster
 ***************************************************************************/
"""

# from osgeo import gdal, ogr, osr
import logging
import rasterio
import numpy as np
import time

from raster_processing import sum_sound_level_3D, resample, reproject, vectorize, check_projection, validate_topo, \
     validate_source_format

formatter = '%(levelname)s: %(asctime)s - %(name)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=formatter, filename='logfile.log') #stream=sys.stdout if stdout, not stderr
logger = logging.getLogger("noise_raster_v1")

### Data paths
src_pths = [r'./data/05111000_SCS_RAST_DAY.ASC', r'./data/05111000_STR_RAST_NGT.ASC']
out_pths = [r'./data/out_test_1.tif', r'./data/out_test_2.tif']

logger.info('First file path: {}'.format(str(src_pths[0])))

### Open source rasters - by using "with open" you don't need to close the files anymore
with rasterio.open(src_pths[0]) as ras1, rasterio.open(src_pths[1]) as ras2:
    logger.info('First file no. of bands: {}'.format(ras1.count))
    logger.info('Second file no. of bands: {}'.format(ras2.count))

    ### Read rasters as numpy 2d arrays
    ras1band = ras1.read(1)
    ras2band = ras2.read(1)
    ### Get affine transformation coefficients from source raster
    out_profile = ras1.profile.copy()

### Read number of rows and cols
numOfRows = ras2band.shape[0] # why do you switch the order 1,2 here for the rows?
numOfRows1 = ras1band.shape[0]
numOfColumns = ras2band.shape[1] # why do you switch the order 1,2 here for the cols?
numOfColumns1 = ras1band.shape[1]

### Print raster information to console
logger.info('Number of Rows in raster array number one: {}'.format(numOfRows1))
logger.info('Number of Rows in raster array number two: {}'.format(numOfRows))
logger.info('Number of Columns in raster array number one: {}'.format(numOfColumns1))
logger.info('Number of Columns in raster array number two: {}'.format(numOfColumns))

### Create input for calc function (3D array with all sound levels stack)
input_3D = np.stack([ras1band, ras2band], axis=0)
# input_3D_test = (1,2,3) - test that also errors are written into the logfile

### Call calculation function
logger.info('Start adding sound levels...')
start = time.time()
data_out = sum_sound_level_3D(input_3D)
end = time.time()
logger.info('Elapsed time [sec] for calculation: {}'.format(end - start))

### Print out information about the returned raster
logger.info('Output array returned from noise level addition: {}'.format(data_out))


###  FROM HERE: we need to discuss if using rasterio is ok or we should go back to gdal ###

# ### Get columns and rows of output array
# cols = data_out.shape[1]
# rows = data_out.shape[0]
#
# ### Create the output raster
# out_affine = out_profile["transform"]
#
# with rasterio.open(out_pths[0], 'w', format='GTiff', affine=out_affine) as dsOut:
#     dsOut.write(data_out)
#
# # driver = gdal.GetDriverByName("GTiff")
# # #dsOut = driver.Create("C:/Users/KL/development/test25832.tif",cols,rows,1,eType=gdal.GDT_Float64)
# # dsOut = driver.Create(out_pth+"/test25832.tif", cols, rows, 1, eType=gdal.GDT_Float64)
# #
# # ### Get affine transformation coefficients from source raster
# # dsOut.SetGeoTransform(ras1.GetGeoTransform())
# #
# # ### Write the output raster
# # bandOut = dsOut.GetRasterBand(1).WriteArray(data_out)
# #
# # ### Set the crs of output raster
# # outRasterSRS = osr.SpatialReference()
# # outRasterSRS.ImportFromEPSG(25832)
# # dsOut.SetProjection(outRasterSRS.ExportToWkt())

