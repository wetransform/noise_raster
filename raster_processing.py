# -*- coding: utf-8 -*-
"""
/***************************************************************************
 NoiseRaster functions
 ***************************************************************************/
"""

from osgeo import gdal, ogr, osr
import os
import numpy as np


def sum_sound_level_3D(sound_levels: np.array):
    """
    INPUT: array of dimension (m,n,l) - stack of sound levels on a 2D map
    OUTPUT: array of dimension (m,n) - final values of cumulative sound levels on a 2D map

    m: y-nodes on the grid
    n: x-nodes on the grid
    l: no. of sound levels per cell's node

    """

    if not isinstance(sound_levels, np.ndarray):
        raise TypeError('Input is not an array')

    if len(sound_levels.shape) != 3:
        raise ValueError('Input array not 3D ')

    l, m, n = sound_levels.shape

    sound_pressures = np.zeros((l, m, n))
    sum_pressures = np.zeros((l, m, n))
    out = np.zeros((l, m, n))

    sound_pressures = np.power(10, 0.1 * sound_levels)
    sum_pressures = np.sum(sound_pressures, axis=0)
    out = 10 * np.log10(sum_pressures)

    return out

def vectorize(out_ds, out_poly):

    #Open raster
    check_ds = gdal.Open(out_ds)

    #Get raster band
    band1 = check_ds.GetRasterBand(1)

    #Set destination SRS
    dest_srs = osr.SpatialReference()
    dest_srs.ImportFromEPSG(3035)

    #Create shapefile
    drv = ogr.GetDriverByName("ESRI Shapefile")
    #Delete existing shapefile if file already exists
    if os.path.exists(out_poly):
        drv.DeleteDataSource(out_poly)
    dst_ds = drv.CreateDataSource(out_poly)
    dst_layer = dst_ds.CreateLayer(out_poly, dest_srs, geom_type=ogr.wkbPolygon)

    #Add field to shapefile
    f1 = ogr.FieldDefn('DN', ogr.OFTInteger)
    dst_layer.CreateField(f1)

    # Vectorize raster
    gdal.Polygonize(srcBand=band1,maskBand=None,outLayer=dst_layer,iPixValField=0)

    # Close dataset
    dst_ds = None
    gdal.Unlink(out_ds)


def check_projection():



def reproject(XRes, YRes, out_ds, ds):

    #Set destination SRS as 3035
    dest_srs = osr.SpatialReference()
    dest_srs.ImportFromEPSG(3035)

    #Set options for reprojection
    wo = gdal.WarpOptions(dstSRS='EPSG:3035', XRes, YRes)

    #Reproject raster
    w = gdal.Warp(destNameOrDestDS=out_ds, srcDSOrSrcDSTab=ds, options=wo)
    
    #Close raster
    w = None

def resample():



def validate_source_format():


def validate_topo():
