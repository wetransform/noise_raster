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
import tempfile
import subprocess
import numpy as np
import logging
import time

# Get current temp directory

cur_temp_dir = tempfile.gettempdir()
rep_temp_dir = cur_temp_dir.replace(os.sep, '/')
temp_dir = rep_temp_dir + "/"

# Initiate logging

logfile = temp_dir + "logfile.txt"
formatter = '%(levelname)s: %(asctime)s - %(name)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=formatter, filename=logfile) #stream=sys.stdout if stdout, not stderr
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


def source_raster_list(*folderpaths):
    """
    INPUT: paths to the raster files
    OUTPUT: list of list of the raster files inside the folders
    """
    rasterlist = []
    out_rasterlist = []

    for path in folderpaths:
        rasterlist.append(path)

    for path in rasterlist:
        rasterfile = [os.path.join(path, f) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

        out_rasterlist.append(rasterfile)
        print(out_rasterlist)
    return out_rasterlist

def check_extent(extent_list:list):
    """
    Check extent of each raster is divisible by 100
    """
    for ds in extent_list:
        # Open dataset
        check_extent = gdal.Open(ds)

        # Get extent of raster
        xmin, xpixel, _, ymax, _, ypixel = check_extent.GetGeoTransform()
        xmin_center = xmin + (xpixel / 2)
        ymax_center = ymax + (abs(ypixel) / 2)

        if (abs(xmin_center % 100) < 0.001) or (abs(ymax_center % 100) < 0.001):
            # Get source file name
            f_name = os.path.basename(ds).split(".")[0]
            # Write file name and extent to logfile
            logger.info('CHECK EXTENT: Filename: {}'.format(str(f_name)))
            logger.info('Pixel value: {}'.format(str(xmin_center)))
            logger.info('Pixel value: {}'.format(str(ymax_center)))
        else:
            pass

def build_virtual_raster(in_vrt:list):
    """
    Build virtual multi-band raster as input to addition
    """
    merged_vrt = temp_dir + 'allnoise.vrt'

    # Set options
    vrto = gdal.BuildVRTOptions(separate=True)
    gdal.BuildVRT(merged_vrt, in_vrt, options=vrto)

    return merged_vrt

    # Write to disk
    merged_vrt = None

def create_zero_array(in_ds):
    """
    Create array as input to energetic addition.
    """
    # Open dataset
    ds = gdal.Open(in_ds)

    data = ds.ReadAsArray()

    # Create masked array with no data values
    zeroData = np.where(data < 0, 0, data)

    return zeroData

def create_raster(sound_array:np.ndarray, out_pth, merged_vrt):
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
    # Delete existing raster if file already exists
    if os.path.exists(out_pth):
        driver.DeleteDataSource(out_pth)
    dsOut = driver.Create(out_pth, cols, rows, 1, eType=gdal.GDT_Float32)

    # Set affine transformation coefficients from source raster
    dsOut.SetGeoTransform(geotransform)

    # Write the output raster
    dsOut.GetRasterBand(1).WriteArray(sound_array)
    dsOut.GetRasterBand(1).SetNoDataValue(0)

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
    elif selectedTableIndex == 1:
        # Lnight
        selectedTable = [49.5, 54.49999, 50, 54.5, 59.49999, 55, 59.5, 64.49999, 60, 64.5, 69.49999, 65, 69.5, '', 70]
    else:
        # Write to logfile
        logger.info('RECLASSIFICATION TABLE WAS NOT SELECTED')
        raise ValueError('Reclassification table not selected')

    # Reclassification output path
    out_reclass = temp_dir + 'reclass.tif'

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


def check_projection(in_data: list):
    """
    Check for existence of a defined projection in GTiff files
    """

    for ds in in_data:

        # Check if raster is asc
        # Convert to lowercase
        src = ds.lower()
        if src.endswith(".asc"):
            # Check is passed for asc files because asc files do not contain crs information.
            pass

        # Check crs definition if data is GTiff
        else:

            # Open dataset
            check_ds = gdal.Open(ds)

            # Check projection of source raster
            prj = check_ds.GetProjection()

            if len(prj) == 0:
                # Get source file name
                f_name = os.path.basename(ds).split(".")[0]
                # Write file name to logfile
                logger.info('SPATIAL REFERENCE NOT DEFINED: Filename: {}'.format(str(f_name)))
                raise ValueError('Spatial reference system is not defined')

def reproject(input_files_path:list):
    """
    Reproject rasters to EPSG:3035. Set no data value to -99.0.
    Source crs of asc files is assumed to be EPSG:25832.
    Source crs of GTiff files is read from the data.
    """

    # List to hold list of reprojected rasters
    reprojectedlist = []

    for input in input_files_path:

        # Check if raster is asc
        # Convert to lowercase
        src = input.lower()
        if src.endswith(".asc"):
            # Warp asc files to EPSG:3035 using EPSG:25832 as defined source crs

            # Get source file name
            f_name = os.path.basename(input).split(".")[0]

            # Create reprojected raster
            out_tif = temp_dir + f_name + "_3035.tif"

            # Reproject rasters to EPSG:3035 in GTiff format, set data type
            gdal.Warp(destNameOrDestDS=out_tif, srcDSOrSrcDSTab=input,
                    options=gdal.WarpOptions(format='GTiff', srcSRS='EPSG:25832', dstSRS='EPSG:3035',
                                            outputType=gdal.GDT_Float32))

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

            # Add raster list to list of lists
            reprojectedlist.append(out_tif)

            # Close raster
            out_tif = None

        else:
            # File is GTiff
            # Warp GTiff files to EPSG:3035 using dataset's defined crs as source crs

            # Get source file name
            f_name = os.path.basename(input).split(".")[0]

            # Create reprojected raster
            out_tif = temp_dir + f_name + "_3035.tif"

            # Open dataset
            check_ds = gdal.Open(input)

            # Check projection of source raster
            prj = check_ds.GetProjection()
            src_srs = osr.SpatialReference(wkt=prj)

            # Reproject rasters to EPSG:3035 in GTiff format, set data type
            gdal.Warp(destNameOrDestDS=out_tif, srcDSOrSrcDSTab=input,
                    options=gdal.WarpOptions(format='GTiff', srcSRS=src_srs, dstSRS='EPSG:3035',
                                            outputType=gdal.GDT_Float32))

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
            out_vrt = temp_dir + str(counter) + "_3035.vrt"

            # Create converted tif
            out_tif = temp_dir + str(counter) + "_3035.tif"

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
        out_vrt = temp_dir + "_3035.vrt"

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
    """

    # Check that the file extension of all input rasters is .asc or .tif
    for srcDs in srcData:
        # Convert to lowercase
        src = srcDs.lower()
        if (src.endswith(".asc")) or (src.endswith(".tif")):
            pass
        else:
            # Get source file name
            f_name = os.path.basename(srcDs).split(".")[0]
            # Write file name to logfile
            logger.info('EXTENSION NOT ASC OR TIF: Filename: {}'.format(str(f_name)))
            raise ValueError('File extension is not .asc or .tif')

def set_nodata_value(in_ds):
    """
    Check the existing no data value of raster.
    Set no data value to -99.0
    """

    # Read existing no data value(s)
    dst_ds = gdal.Open(in_ds, gdal.GA_Update)
    ndv = dst_ds.GetRasterBand(1).GetNoDataValue()
    newndv = -99.0
    band1 = dst_ds.GetRasterBand(1).ReadAsArray()
    band1[band1 == ndv] = newndv

    # Set no data value for source rasters to -99.0
    dst_ds.GetRasterBand(1).SetNoDataValue(newndv)
    dst_ds.GetRasterBand(1).WriteArray(band1)
    dst_ds = None

    return in_ds