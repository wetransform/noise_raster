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

def vectorize(in_ds, out_poly):

    #Open raster
    check_ds = gdal.Open(in_ds)

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
    check_ds = None
    dst_ds = None
    gdal.Unlink(out_ds)


def check_projection(in_ds):

    ###Open raster
    check_ds = gdal.Open(in_ds)

    ###Check projection of source raster
    prj = check_ds.GetProjection()
    src_srs = osr.SpatialReference(wkt=prj)
    epsg_code = src_srs.GetAttrValue('AUTHORITY', 1)
    print(epsg_code)
    
    if src_srs is None:
        raise ValueError('Spatial reference system is not defined')



def reproject(XRes, YRes, out_ds, ds):

    #Set destination SRS as 3035
    dest_srs = osr.SpatialReference()
    dest_srs.ImportFromEPSG(3035)

    #Set options for reprojection
    wo = gdal.WarpOptions(dstSRS='EPSG:3035', xRes=XRes, yRes=YRes)

    #Reproject raster
    w = gdal.Warp(destNameOrDestDS=out_ds, srcDSOrSrcDSTab=ds, options=wo)

    #Close raster
    w = None

def merge_rasters(ras_list:list, out_vrt):

    ###Check resolution of input rasters
    #for ras in ras_list:
        ###Open raster
        #r = gdal.Open(ras)
        #gt = r.GetGeoTransform()
        ###Get cell resolution
        #pixelSizeX = gt[1]
        #pixelSizeY = -gt[5]
        ###Check size resolution is 10
        #if not pixelSizeX == pixelSizeY != 10:
            #raise ValueError('Cell resolution is not 10')

    ###Merge rasters in list
    gdal.BuildVRT(out_vrt, ras_list)
    out_vrt = None

def validate_source_format():
    pass

def validate_topo():
    pass