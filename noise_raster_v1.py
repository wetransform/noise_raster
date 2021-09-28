# -*- coding: utf-8 -*-
"""
/***************************************************************************
 NoiseRaster
 ***************************************************************************/
"""

from osgeo import gdal, ogr, osr
import os
import numpy as np

from raster_processing import sum_sound_level_3D, resample, reproject, vectorize, check_projection, validate_topo, validate_source_format

###Get source data paths
src_pths = [r"C:\Users\KL\Downloads\Beispieldaten_DF48_Raster\Beispieldaten_DF48_Raster\Rail_NW\Daten_UBA\Daten_UBA\Raster_DAY_SCS\DAY_SCS_RAST_DAY.ASC",r"C:\Users\KL\Downloads\Beispieldaten_DF48_Raster\Beispieldaten_DF48_Raster\Road_NW\Projekt_ruhige_Gebiete_Rasterlärmkarten_alle_Strassen\Projekt_ruhige_Gebiete_Rasterlärmkarten_alle_Strassen\05111000_STR_RAST_NGT.ASC\DAY_STR_RAST_NGT.ASC"]

###Open source rasters
ras1 = gdal.Open(src_pths[0])
ras2 = gdal.Open(src_pths[1])
print(str(src_pths[0]))

##Read rasters as numpy 2d arrays
ras1band = np.array(ras1.GetRasterBand(1).ReadAsArray())
ras2band = np.array(ras2.GetRasterBand(1).ReadAsArray())

###Read number of rows
numOfRows = ras2band.shape[0]
numOfRows1 = ras1band.shape[0]

##Print raster information to console
print('Number of Rows in raster array number one : ', numOfRows1)
print('Number of Rows in raster array number two : ', numOfRows)

###Read number of columns
numOfColumns = ras2band.shape[1]
numOfColumns1 = ras1band.shape[1]

##Print raster information to console
print('Number of Columns in raster array number one: ', numOfColumns1)
print('Number of Columns in raster array number two : ', numOfColumns)

### Create input for calc function (3D array with all sound levels stack)
input_3D = np.stack([ras1band, ras2band], axis=0)

###Call calculation function

data_out = sum_sound_level_3D(input_3D)

###Print out information about the returned raster
print('Array data_out returned from calculation:', data_out)

###Get columns and rows of output array
cols = data_out.shape[1]
rows = data_out.shape[0]

####Create the output raster
driver = gdal.GetDriverByName("GTiff")
dsOut = driver.Create("C:/Users/KL/development/test25832.tif",cols,rows,1,eType=gdal.GDT_Float64)

###Get affine transformation coefficients from source raster
dsOut.SetGeoTransform(ras1.GetGeoTransform())

###Write the output raster
bandOut = dsOut.GetRasterBand(1).WriteArray(data_out)

###Set the crs of output raster
outRasterSRS = osr.SpatialReference()
outRasterSRS.ImportFromEPSG(25832)
dsOut.SetProjection(outRasterSRS.ExportToWkt())

###Close the datasets
ras1 = None
ras2 = None
ras1band = None
ras2band = None
bandOut = None
dsOut = None
