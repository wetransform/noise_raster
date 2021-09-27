# -*- coding: utf-8 -*-
"""
/***************************************************************************
 NoiseRaster
 ***************************************************************************/
"""
from typing import Any

from osgeo import gdal, ogr, osr
import os
import glob
import numpy as np

from raster_processing import sum_sound_level_3D, merge_rasters, reproject, vectorize, check_projection, validate_source_format

###Input rail noise folder

src_pths_ra =r'C:\Users\KL\Downloads\Beispieldaten_DF48_Raster\Beispieldaten_DF48_Raster\Rail_EBA\RLK_SchVK_AllSrcs - D\*LDEN*.asc'

###Input road noise folder path

src_pths_ro =r'C:\Users\KL\Downloads\Beispieldaten_DF48_Raster\Beispieldaten_DF48_Raster\Road_NW\Projekt_ruhige_Gebiete_Rasterlärmkarten_alle_Strassen\Projekt_ruhige_Gebiete_Rasterlärmkarten_alle_Strassen\*DEN*.ASC'

###Make list of all input rail noise source rasters
src_ra = glob.glob(src_pths_ra)

# src_ra = []
# for root, directories, file in os.walk(src_pths_ra):
#     for file in file:
#         if file.endswith(".asc"):
#             src_ra.append(os.path.join(root,file))

###Make list of all input road noise source rasters
src_ro = glob.glob(src_pths_ro)

# src_ro = []
# for root, directories, file in os.walk(src_pths_ro):
#     for file in file:
#         if file.endswith(".ASC"):
#             src_ro.append(os.path.join(root,file))

###Check CRS of each input raster

#for ras_ra in src_ra:
    #check_projection(ras_ra)

#for ras_ro in src_ro:
    #check_projection(ras_ro)

###Set no data value of each input raster

#validate_source_format(src_ra)

#validate_source_format(src_ro)

###Create output tif raster for merged rail noise

#out_tif_ra = "C:/Users/KL/development/test_ra.tif"

###Create output tif raster for merged road noise

out_tif_ro = "C:/Users/KL/development/test_ro_v2.tif"

###Merge all input rail noise rasters

#merge_rasters(src_ra,out_tif_ra)

###Merge all input road noise rasters

merge_rasters(src_ro,out_tif_ro)

###Open merged source rasters
# ras1 = gdal.Open(out_tif_ro)
# ras2 = gdal.Open(out_tif_ra)
#
# ##Read rasters as numpy 2d arrays
# ras1band = np.array(ras1.GetRasterBand(1).ReadAsArray())
# ras2band = np.array(ras2.GetRasterBand(1).ReadAsArray())
#
# ###Read number of rows
# numOfRows = ras2band.shape[0]
# numOfRows1 = ras1band.shape[0]
#
# ##Print raster information to console
# print('Number of Rows in raster array number one : ', numOfRows1)
# print('Number of Rows in raster array number two : ', numOfRows)
#
# ###Read number of columns
# numOfColumns = ras2band.shape[1]
# numOfColumns1 = ras1band.shape[1]
#
# ##Print raster information to console
# print('Number of Columns in raster array number one: ', numOfColumns1)
# print('Number of Columns in raster array number two : ', numOfColumns)

### Create input for calc function (3D array with all sound levels stack)
###input_3D = np.stack([ras1band, ras2band], axis=0)

###Call calculation function

###data_out = sum_sound_level_3D(input_3D)

###Print out information about the returned raster
#print('Array data_out returned from calculation:', data_out)

###Get columns and rows of output array
# cols = data_out.shape[1]
# rows = data_out.shape[0]
#
# ####Create the output raster
# driver = gdal.GetDriverByName("GTiff")
# dsOut = driver.Create("C:/Users/KL/development/test3035.tif",cols,rows,1,eType=gdal.GDT_Float64)
#
# ###Get affine transformation coefficients from source raster
# dsOut.SetGeoTransform(ras1.GetGeoTransform())
#
# ###Write the output raster
# bandOut = dsOut.GetRasterBand(1).WriteArray(data_out)
#
# ###Set the crs of output raster
# outRasterSRS = osr.SpatialReference()
# outRasterSRS.ImportFromEPSG(3035)
# dsOut.SetProjection(outRasterSRS.ExportToWkt())

###Close the datasets
#bandOut = None
#dsOut = None

###Vectorize energetically added raster including all noise sources
#out_poly_ra = "C:/Users/KL/development/test_ra_v2.shp"
out_poly_ro = "C:/Users/KL/development/test_ro_v2.shp"
#vectorize(out_tif_ra, out_poly_ra)
vectorize(out_tif_ro, out_poly_ro)