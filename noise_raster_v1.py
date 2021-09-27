# -*- coding: utf-8 -*-
"""
/***************************************************************************
 NoiseRaster
 ***************************************************************************/
"""
from qgis.PyQt.QtCore import QSettings, QTranslator, QCoreApplication
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QAction
from qgis.core import QgsProject
from osgeo import gdal, ogr, osr
import os
import numpy as np

    ###get source data paths
    ###src_pths = self.dlg.mQgsFileWidget.splitFilePaths(self.dlg.mQgsFileWidget.filePath())
    src_pths = ["/05111000_SCS_RAST_DAY.ASC", "/05111000_STR_RAST_NGT.ASC"]
    ###Open source rasters
    ras1 = gdal.Open(src_pths[0])
    ras2 = gdal.Open(src_pths[1])
	##Read rasters as numpy 2d arrays
    ras1band = np.array(ras1.GetRasterBand(1).ReadAsArray())
    ras2band = np.array(ras2.GetRasterBand(1).ReadAsArray())
    numOfRows = ras2band.shape[0]
    numOfRows1 = ras1band.shape[0]
	##Print raster information to console
    print('Number of Rows in raster array number one : ', numOfRows1)
    print('Number of Rows in raster array number two : ', numOfRows)
    numOfColumns = ras2band.shape[1]
    numOfColumns1 = ras1band.shape[1]
    print('Number of Columns in raster array number one: ', numOfColumns1)
    print('Number of Columns in raster array number two : ', numOfColumns)
    #calculation
    grid = []
    input_list = [ras1band, ras2band]
    for arr in input_list:
        sum_sound_levels = 10*np.log10(np.sum(np.power(10, 0.1*arr)))
        grid.append(sum_sound_levels)
    data_out = np.array(grid)
    print('Array data_out returned from calculation:', data_out)
    numOfRows11 = data_out.shape[0]
    numOfColumns11 = data_out.shape[1]
    print('Number of Rows : ', numOfRows11)
    print('Number of Columns : ', numOfColumns11)
    # Write the out file
    driver = gdal.GetDriverByName("GTiff")
    dsOut = driver.Create("C:/Users/KL/development/test3035.tif",1770,2570,1,eType=gdal.GDT_Float64)
    bandOut = dsOut.GetRasterBand(1).WriteArray(data_out)
    # Close the datasets
    ras1 = None
    ras2 = None
    ras1band = None
    ras2band = None
    bandOut = None
    dsOut = None
