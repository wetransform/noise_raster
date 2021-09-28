# -*- coding: utf-8 -*-
"""
/***************************************************************************
 NoiseRaster
 ***************************************************************************/
"""

from osgeo import gdal, ogr, osr
from sound_level_add_3D_efficient import sum_sound_level_3D 
import os
import numpy as np

from .raster_processing import sum_sound_level_3D, resample, reproject, vectorize, check_projection, validate_topo, validate_source_format

    ###Get source data paths
    src_pths = ["/05111000_SCS_RAST_DAY.ASC", "/05111000_STR_RAST_NGT.ASC"]

    ###Open source rasters
    ras1 = gdal.Open(src_pths[0])
    ras2 = gdal.Open(src_pths[1])

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

	
    ####Calculation - you don't need this anymore
    ## grid = []
    ## input_list = [ras1band, ras2band]

    ### Create input for calc function (3D array with all sound levels stack) 
    input_3D = np.stack([ras1band, ras2band], axis=0)

    ###Call calculation function
    ## sum_sound_level_3D(input_list)
    out = sum_sound_level_3D(input_3D) 
    ## check if out is what you expect to be	


    ###Print out information about the returned raster
    print('Array data_out returned from calculation:', data_out)
    numOfRows11 = data_out.shape[0]
    numOfColumns11 = data_out.shape[1]
    print('Number of Rows : ', numOfRows11)
    print('Number of Columns : ', numOfColumns11)

    ####Write the output raster
    driver = gdal.GetDriverByName("GTiff")
    dsOut = driver.Create("C:/Users/KL/development/test3035.tif",1770,2570,1,eType=gdal.GDT_Float64)
    bandOut = dsOut.GetRasterBand(1).WriteArray(data_out)

    ###Close the datasets
    ras1 = None
    ras2 = None
    ras1band = None
    ras2band = None
    bandOut = None
    dsOut = None
