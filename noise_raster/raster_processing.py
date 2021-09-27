# -*- coding: utf-8 -*-
"""
/***************************************************************************
 NoiseRaster functions
 ***************************************************************************/
"""

from osgeo import gdal, ogr, osr
from qgis.core import QgsRasterLayer, QgsProcessing
import sys
import os
import processing
import subprocess
import numpy as np
import logging
import time

# Initiate logging

formatter = '%(levelname)s: %(asctime)s - %(name)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=formatter, filename='C:/Users/KL/development/logfile.txt') #stream=sys.stdout if stdout, not stderr
logger = logging.getLogger("noise_raster_v1")


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

def source_raster_list(folderpath, folderpath2=None, folderpath3=None):
    """
    Create list of list of rasters in each source folder
    """
    rasterlist = [folderpath]
    out_rasterlist = []

    if len(folderpath2) != 0:
        rasterlist.append(folderpath2)

    if len(folderpath3) != 0:
        rasterlist.append(folderpath3)

    for l in rasterlist:
        newl = [os.path.join(l, f) for f in os.listdir(l) if os.path.isfile(os.path.join(l, f))]

        out_rasterlist.append(newl)

    return out_rasterlist

def check_extent(extent_list:list):
    """
    Check extent of each raster is divisible by 100
    """
    for ds in extent_list:
        # Open dataset
        check_extent = gdal.Open(ds)

        xmin, xpixel, _, ymax, _, ypixel = check_extent.GetGeoTransform()
        width, height = check_extent.RasterXSize, check_extent.RasterYSize
        xmax = xmin + width
        ymin = ymax + height

def build_virtual_raster(in_vrt:list):
    """
    Build virtual multi-band raster as input to addition
    """
    merged_vrt = 'C:/Users/KL/development/all_noise_vrt.vrt'

    # Set options
    vrto = gdal.BuildVRTOptions(separate=True)
    gdal.BuildVRT(merged_vrt, in_vrt, options=vrto)

    return merged_vrt

    # Write to disk
    merged_vrt = None

def create_masked_array(in_ds):
    """
    Create masked array as input to energetic addition. Masked arrays handle no data values.
    """
    # Open dataset
    ds = gdal.Open(in_ds)

    data = ds.ReadAsArray()

    # Create masked array with no data values
    maskedData = np.ma.array(data, mask=(data == -99.0))

    return maskedData

def create_raster(sound_array, out_pth, merged_vrt):
    """
    Create raster based on energetically added array
    """

    # Get geotranform of merged_vrt
    ds = gdal.Open(merged_vrt)
    geotransform = ds.GetGeoTransform()

    # Get columns and rows of output array
    cols = sound_array.shape[1]
    rows = sound_array.shape[0]

    # Create the output raster
    driver = gdal.GetDriverByName("GTiff")
    dsOut = driver.Create(out_pth, cols, rows, 1, eType=gdal.GDT_Float32)

    # Set affine transformation coefficients from source raster
    dsOut.SetGeoTransform(geotransform)

    # Write the output raster
    dsOut.GetRasterBand(1).WriteArray(sound_array)
    dsOut.GetRasterBand(1).SetNoDataValue(-99.0)

    # Set the crs of output raster
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(3035)
    dsOut.SetProjection(outRasterSRS.ExportToWkt())

    return out_pth

    # Close the datasets
    dsOut = None
    out_pth = None

def vectorize(in_ds, out_poly, selectedTableIndex):
    """
    Reclassify raster.
    Convert reclassified raster to polygon.
    """

    # Set destination SRS of target shapefile
    dest_srs = osr.SpatialReference()
    dest_srs.ImportFromEPSG(3035)

    # Create shapefile
    drv = ogr.GetDriverByName("ESRI Shapefile")
    # Delete existing shapefile if file already exists
    if os.path.exists(out_poly):
        drv.DeleteDataSource(out_poly)
    dst_ds = drv.CreateDataSource(out_poly)
    dst_layer = dst_ds.CreateLayer(out_poly, dest_srs, geom_type=ogr.wkbPolygon)

    # Add field to shapefile
    noise_field = ogr.FieldDefn('Noise', ogr.OFTReal)
    dst_layer.CreateField(noise_field)

    # Apply selected reclassification values from combobox
    if selectedTableIndex == 0:
        # Lden
        selectedTable = [54.5, 59.49999, 55, 59.5, 64.49999, 60, 64.5, 69.49999, 65, 69.5, 74.49999, 70, 74.5, '', 75]
    else:
        # Lnight
        selectedTable = [49.5, 54.49999, 50, 54.5, 59.49999, 55, 59.5, 64.49999, 60, 64.5, 69.49999, 65, 69.5, '', 70]

    # Reclassification output path
    out_reclass = 'C:/Users/KL/development/reclass_vrt.tif'

    # Set reclassify algorithm parameters
    alg_params = {
        'DATA_TYPE': 5,
        'INPUT_RASTER': in_ds,
        'NODATA_FOR_MISSING': False,
        'NO_DATA': -99.0,
        'RANGE_BOUNDARIES': 2,
        'RASTER_BAND': 1,
        'TABLE': selectedTable,
        'OUTPUT': out_reclass
    }

    # Reclassify raster
    result = processing.run('native:reclassifybytable', alg_params)

    # Open reclassified raster
    check_ds = gdal.Open(result['OUTPUT'])

    # Get reclassified raster band
    band1 = check_ds.GetRasterBand(1)

    # Vectorize reclassified raster to create polygon noise contours
    gdal.Polygonize(srcBand=band1, maskBand=None, outLayer=dst_layer, iPixValField=0)

    # Close dataset
    check_ds = None
    dst_ds = None


def check_projection(in_data:list):
    """
    Check for existence of a defined projection
    """

    # Open raster
    for ds in in_data:
        check_ds = gdal.Open(ds)

        # Check projection of source raster
        prj = check_ds.GetProjection()
        src_srs = osr.SpatialReference(wkt=prj)
        epsg_code = src_srs.GetAttrValue('AUTHORITY', 1)

        if src_srs is None:
            raise ValueError('Spatial reference system is not defined')

def reproject(input_files_path:list):
    """
    Reproject rasters to EPSG:3035
    """

    # List to hold list of reprojected rasters
    reprojectedlist = []

    for input in input_files_path:
        # Get source file name
        f_name = os.path.basename(input).split(".")[0]

        # Create reprojected raster
        out_tif = "C:/Users/KL/development/" + f_name + "_3035.tif"

        # Reproject rasters
        gdal.Warp(destNameOrDestDS=out_tif, srcDSOrSrcDSTab=input, options=gdal.WarpOptions(format='GTiff', dstSRS='EPSG:3035', outputType=gdal.GDT_Float32, srcNodata=-99.0, dstNodata=-99.0))

        # Add raster list to list of lists
        reprojectedlist.append(out_tif)

        # Close raster
        out_tif = None

    return reprojectedlist

def merge_rasters(input_files_path:list, out=None):
    """
    Merge each list of rasters into a single raster
    """

    # Create list of merged rasters
    merged_ras = []

    if len(input_files_path) != 0:
        counter = 1

        for input in input_files_path:

            # Create merged vrt
            out_vrt = "C:/Users/KL/development/" + str(counter) + "_3035.vrt"

            # Create converted tif
            out_tif = "C:/Users/KL/development/" + str(counter) + "_3035.tif"

            # Set vrt options
            gdal.BuildVRT(out_vrt, input, resolution='highest', resampleAlg=gdal.gdalconst.GRA_Max, outputSRS='EPSG:3035', srcNodata=-99.0)

            # Set translate options
            to = gdal.TranslateOptions(format="GTiff", outputSRS="EPSG:3035", noData=-99.0, outputType=gdal.GDT_Float32)

            # Convert vrt file to tif
            gdal.Translate(out_tif, out_vrt, options=to)

            # Add merged raster to list of merged rasters
            merged_ras.append(out_vrt)

            # Close raster
            out_vrt = None

            # Add 1 to counter
            counter += 1

    else:
        """
        Merge only one list of rasters
        """

        # Create merged vrt
        out_vrt = "C:/Users/KL/development/vrt_3035.vrt"

        # Set vrt options
        gdal.BuildVRT(out_vrt, input, resolution='highest', resampleAlg=gdal.gdalconst.GRA_Max, outputSRS='EPSG:3035', srcNodata=-99.0)

        # Set translate options
        to = gdal.TranslateOptions(format="GTiff", outputSRS="EPSG:3035", noData=-99.0, outputType=gdal.GDT_Float32)

        # Convert vrt file to tif
        gdal.Translate(out, out_vrt, options=to)

        # Close raster
        out_vrt = None

    return merged_ras

def validate_source_format(srcData:list):
    """
    Check that the file extension of all input rasters is .asc or .tif.
    Ensure all input files have the same no data value and data type.
    """

    # Check that the file extension of all input rasters is .asc or .tif
    for srcDs in srcData:
        # Convert to lowercase
        src = srcDs.lower()
        if (src.endswith(".asc")) or (src.endswith(".tif")):
            pass
        else:
            raise ValueError('File extension is not .asc or .tif')

    # Create list to hold lists of translated rasters
    out_ras = []

    for srcDs in srcData:
        # Ensure all input files have the same no data value and data type

        # Get source file name
        f_name = os.path.basename(srcDs).split(".")[0]

        # Read all negative pixel values
        # Open raster
        n_pixels = gdal.Open(srcDs)

        # Read raster
        band = n_pixels.GetRasterBand(1).ReadAsArray()

        # Loop pixel values in array and write negative values to file
        for arr in band:
            for pixel in arr:
                if 0 > pixel > -90:
                    # Start logging raster information to file
                    logger.info('Raster dataset: {}'.format(str(srcDs)))
                    logger.info('Pixel value: {}'.format(str(pixel)))

        # Set all negative pixel values to 0. Recommendation from https://gdal.org/programs/gdal_calc.html
        out_calc = "C:/Users/KL/development/" + f_name + "calc.tif"

        # Parameters for gdal raster calculator
        alg_params = {
            'BAND_A': 1,
            'FORMULA': 'A*(A>0)',
            'INPUT_A': srcDs,
            'NO_DATA': 0,
            'RTYPE': 5,
            'OUTPUT': out_calc
                    }
        # Execute gdal raster calculator
        result= processing.run('gdal:rastercalculator', alg_params)

        # Set translate options
        to = gdal.TranslateOptions(format="GTiff", outputType=gdal.GDT_Float32)

        # Convert files to tif with same data type and same no data value
        out_tif = "C:/Users/KL/development/" + f_name + ".tif"
        gdal.Translate(out_tif, result['OUTPUT'], options=to)

        # Read existing no data value(s)
        dst_ds = gdal.Open(out_tif, gdal.GA_Update)
        ndv = dst_ds.GetRasterBand(1).GetNoDataValue()
        newndv = -99.0
        band1 = dst_ds.GetRasterBand(1).ReadAsArray()
        band1[band1 == ndv] = newndv

        # Set no data value for source rasters to -99.0
        dst_ds.GetRasterBand(1).SetNoDataValue(newndv)
        dst_ds.GetRasterBand(1).WriteArray(band1)
        dst_ds = None

        # Add tifs with unified no data value to list
        out_ras.append(out_tif)

    return out_ras