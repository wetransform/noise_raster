# -*- coding: utf-8 -*-
"""
/***************************************************************************
 NoiseRaster functions
 ***************************************************************************/
"""

from osgeo import gdal, ogr, osr
from qgis.core import QgsRasterLayer, QgsProcessing
import os
import shutil
import processing
import tempfile
import numpy as np
import logging
from datetime import datetime
import noise_raster.constants as c

# Reclassification tables
"""
The selected table values are required input for the "Reclassify by Table" operation. The user selects "Lden" or "Lnight" in the user interface of the tool. 
"Lden" categories are used for daytime noise and "Lnight" are used for night time noise. These value ranges are standards used in noise mapping. 
The input raster is reclassified so that each raster cell gets the same value if it's cell value falls into one of the given ranges. 
This pre-processing step dramatically improves the performance of the vectorization process.
"""
# Lden
selectedTableLden = [54.5, 59.49999, 55, 59.5, 64.49999, 60, 64.5, 69.49999, 65, 69.5, 74.49999, 70, 74.5, '', 75]

# Lnight
selectedTableLnight = [49.5, 54.49999, 50, 54.5, 59.49999, 55, 59.5, 64.49999, 60, 64.5, 69.49999, 65, 69.5, '', 70]

def progress_callback(complete, message, callback_data):
    """
    Emit progress report in numbers for 10% intervals and dots for 3%
    https://stackoverflow.com/questions/68025043/adding-a-progress-bar-to-gdal-translate
    """
    if int(complete*100) % 10 == 0:
        print(f'{complete*100:.0f}', end='', flush=True)
    elif int(complete*100) % 3 == 0:
        print(f'{callback_data}', end='', flush=True)

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

    sound_levels.astype(np.float32)

    l, m, n = sound_levels.shape

    sound_pressures = np.zeros((l, m, n)).astype(np.float32)
    sum_pressures = np.zeros((l, m, n)).astype(np.float32)
    out = np.zeros((l, m, n)).astype(np.float32)

    sound_pressures = np.power(10, 0.1 * sound_levels)
    sum_pressures = np.sum(sound_pressures, axis=0)
    out = (10 * np.log10(sum_pressures)).round(decimals=1)

    return out

def create_temp_directory():
    """
    Create a sub-directory in the current temp folder to hold intermediate files.
    """
    # Get current date and time
    DATETIME = datetime.now()
    # Convert date, time to string in format: dd_mm_YY_H_M_S
    DATETIME_str = DATETIME.strftime("%d_%m_%Y_%H_%M_%S")
    # Create sub directory named 'noise_<current date and time>' in temp folder to collect intermediate files
    os.makedirs(c.root_temp_dir + "/noise_" + DATETIME_str)
    temp_dir = c.root_temp_dir + "/noise_" + DATETIME_str + "/"

    return temp_dir, DATETIME_str

global logger

def start_logging():
    """
    Initiate logging
    """

    logfile = c.root_temp_dir + "/logfile.txt"
    formatter = '%(levelname)s: %(asctime)s - %(name)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=formatter,
                        filename=logfile)  # stream=sys.stdout if stdout, not stderr
    global logger
    logger = logging.getLogger("noise_raster_v1")

    return logger

def log(message, variable=None):
    """
    Write message to log file
    """
    if variable:
        logger.info(message.format(variable))
    else:
        logger.info(message)

def log_console(message, variable=None):
    """
    Write message to python console and log file
    """
    log(message, variable)
    if variable:
        print(message.format(variable), end='', flush=True)
    else:
        print(message, end='', flush=True)


def source_raster_list(*folderpaths):
    """
    INPUT: paths to the raster files
    OUTPUT: list of list of the raster files inside the folders
    """
    rasterlist = []
    out_rasterlist = []

    for path in folderpaths:
        if len(path) != 0:
            rasterlist.append(path)

    for path in rasterlist:
        rasterfile = [os.path.join(path, f) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

        out_rasterlist.append(rasterfile)
        
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
            log('CHECK EXTENT: Filename: {}', (str(f_name)))
            log('Pixel value: {}', (str(xmin_center)))
            log('Pixel value: {}', (str(ymax_center)))
            # Close dataset
            check_extent = None
        else:
            # Close dataset
            check_extent = None
            pass



def build_virtual_raster(in_vrt:list, temp_dir):
    """
    Build virtual multi-band raster as input to addition
    """
    merged_vrt = temp_dir + 'allnoise.vrt'

    # Set options
    vrto = gdal.BuildVRTOptions(separate=True)
    gdal.BuildVRT(merged_vrt, in_vrt, options=vrto)

    return merged_vrt

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

def create_raster(sound_array:np.ndarray, merged_vrt, out_pth=None):
    """
    Create raster based on energetically added array
    """

    # Add tif name and file extension to directory file path.
    out_pth_ext = os.path.join(out_pth, c.REPROJECTED_TIF25832)

    # Get geotranform of merged_vrt
    ds = gdal.Open(merged_vrt)
    geotransform = ds.GetGeoTransform()

    # Get columns and rows of output array
    cols = sound_array.shape[1]
    rows = sound_array.shape[0]

    # Create the output raster
    driver = gdal.GetDriverByName("GTiff")
    dsOut = driver.Create(out_pth_ext, cols, rows, 1, eType=gdal.GDT_Float32)

    # Set affine transformation coefficients from source raster
    dsOut.SetGeoTransform(geotransform)

    # Write the output raster
    dsOut.GetRasterBand(1).WriteArray(sound_array)
    dsOut.GetRasterBand(1).SetNoDataValue(0)

    # Set the crs of output raster
    # TODO Make configurable to support sources with other reference systems
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(25832)
    dsOut.SetProjection(outRasterSRS.ExportToWkt())

    return out_pth_ext

    # Close the datasets
    dsOut = None
    out_pth_ext = None
    ds = None

def vectorize(in_ds, out_poly, selectedTableIndex, temp_dir):
    """
    Create polygon for vectorization. Reclassify source input raster.
    Vectorize reclassified raster.
    Reproject polygon to EPSG:3035.
    """

    # Add shp name and file extension to directory path
    out_poly_pth = os.path.join(out_poly, c.REPROJECTED_SHP25832)

    # Set destination SRS of target shapefile
    dest_srs = osr.SpatialReference()
    dest_srs.ImportFromEPSG(25832)

    # Create shapefile
    drv = ogr.GetDriverByName("ESRI Shapefile")
    # Delete existing shapefile if file already exists
    if os.path.exists(out_poly_pth):
        drv.DeleteDataSource(out_poly_pth)
    dst_ds = drv.CreateDataSource(out_poly_pth)
    dst_layer = dst_ds.CreateLayer(out_poly_pth, dest_srs, geom_type=ogr.wkbPolygon)

    # Add field to shapefile
    noise_field = ogr.FieldDefn('Noise', ogr.OFTReal)
    dst_layer.CreateField(noise_field)

    # Apply selected reclassification values from combobox
    if selectedTableIndex == 0:
        # Lden
        selectedTable = selectedTableLden
    elif selectedTableIndex == 1:
        # Lnight
        selectedTable = selectedTableLnight
    else:
        # Write to logfile
        log('RECLASSIFICATION TABLE WAS NOT SELECTED')
        raise ValueError('Reclassification table not selected')

    # Reclassification output path
    out_reclass = temp_dir + 'reclass.tif'

    # Set reclassify algorithm parameters
    # MISSING values that do not fall in the lden and lnight value ranges are classified as No Data and given a value of 0
    alg_params = {
        'DATA_TYPE': 5,
        'INPUT_RASTER': in_ds,
        'NODATA_FOR_MISSING': True,
        'NO_DATA': 0,
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

    # Write progress to python console
    log_console('\nRunning GDAL polygonize...')

    # Vectorize reclassified raster to create polygon noise contours
    """
    Values that fall outside of the lden and lnight value ranges should not be carried over into the shapefile polygon that is created. 
    In order to remove these values, the no data value is set to 0 for the reclassified raster and then the raster itself is used as the raster mask 
    (instead of creating a separate raster mask file). When calling Polygonize, setting the maskBand parameter to the raster itself 
    removes the no data values identified as 0 from the polygon that is created.
    """
    gdal.Polygonize(srcBand=band1, maskBand=band1, outLayer=dst_layer, iPixValField=0, callback=progress_callback, callback_data='.')

    # Close dataset
    check_ds = None
    dst_ds = None

    # Reproject polygon to EPSG:3035
    out_poly_rprj = os.path.join(out_poly, c.REPROJECTED_SHP3035)

    # Set reprojection parameters

    alg_params_shp = {
        'INPUT': out_poly_pth,
        'TARGET_CRS': 'EPSG:3035',
        'OUTPUT': out_poly_rprj
    }

    # Reproject polygon
    processing.run('native:reprojectlayer', alg_params_shp)

    # Close dataset
    out_poly_rprj = None


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
                log('Filename: {}', (str(f_name)))
                raise ValueError('Spatial reference system is not defined')

            check_ds = None

def reproject(input_files_path:list, temp_dir):
    """
    Reproject tifs to EPSG:25832. Source crs of GTiff files is read from the data.
    Translate asc files to tifs. Source crs of asc files is assumed to be EPSG:25832.
    Set no data value to -99.0.
    """

    # List to hold list of reprojected rasters
    reprojectedlist = []

    for input in input_files_path:

        # Get source directory name
        d_path = os.path.dirname(input)
        d_name = os.path.basename(d_path)

        # Get source file name
        f_name = os.path.basename(input).split(".")[0]

        # Check if raster is asc
        # Convert to lowercase
        src = input.lower()
        if src.endswith(".asc"):
            # Translate asc files to tif using EPSG:25832 as defined source crs

            # Create tif raster
            out_tif = temp_dir + d_name + "_" + f_name + "_25832.tif"

            # Set translate options
            to = gdal.TranslateOptions(format="GTiff", outputSRS='EPSG:25832', outputType=gdal.GDT_Float32)

            # Convert asc file to tif
            gdal.Translate(out_tif, input, options=to)

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
            # Warp GTiff using dataset's defined crs as source crs

            # Create reprojected raster
            out_tif = temp_dir + d_name + "_" + f_name + "_25832.tif"

            # Open dataset
            check_ds = gdal.Open(input)

            # Check projection of source raster
            prj = check_ds.GetProjection()
            src_srs = osr.SpatialReference(wkt=prj)

            # Reproject rasters to EPSG:25832 in GTiff format, set data type
            gdal.Warp(destNameOrDestDS=out_tif, srcDSOrSrcDSTab=input,
                    options=gdal.WarpOptions(format='GTiff', srcSRS=src_srs, dstSRS='EPSG:25832',
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

def merge_rasters(input_files_path:list, temp_dir, out_pth=None):
    """
    Merge lists of rasters into a single raster or merge a single list of rasters to a single raster.
    """

    # Create list of merged rasters
    merged_ras = []

    if len(input_files_path) > 1:
        counter = 1

        for input in input_files_path:

            # Create merged vrt
            out_vrt = temp_dir + str(counter) + "_25832.vrt"

            # Create converted tif
            out_tif = temp_dir + str(counter) + "_25832.tif"

            # Write progress to console
            log_console("\nRunning GDAL BuildVRT to create merged raster based on list of rasters...")

            # Set vrt options
            gdal.BuildVRT(out_vrt, input, resolution='highest', resampleAlg=gdal.gdalconst.GRA_Max, outputSRS='EPSG:25832', srcNodata=-99.0, callback=progress_callback, callback_data='.')

            # Write progress to console
            log_console("\nRunning GDAL Translate to convert virtual merged raster to tif...")

            # Set translate options
            to = gdal.TranslateOptions(format="GTiff", outputSRS="EPSG:25832", noData=-99.0, outputType=gdal.GDT_Float32, callback=progress_callback, callback_data='.')

            # Convert vrt file to tif
            gdal.Translate(out_tif, out_vrt, options=to)

            # Add merged raster to list of merged rasters
            merged_ras.append(out_vrt)

            # Close raster
            out_vrt = None
            out_tif = None

            # Add 1 to counter
            counter += 1

        return merged_ras

    else:
        """
        Merge only one list of rasters
        """

        # Create merged vrt
        out_vrt = temp_dir + "singleList_25832.vrt"

        # Write progress to console
        log_console("\nRunning GDAL BuildVRT to create merged raster based on list of rasters...")

        # Set vrt options
        gdal.BuildVRT(out_vrt, input_files_path[0], resolution='highest', resampleAlg=gdal.gdalconst.GRA_Max, outputSRS='EPSG:25832', srcNodata=-99.0, callback=progress_callback, callback_data='.')

        # Write progress to console
        log_console("\nRunning GDAL Translate to convert virtual merged raster to tif...")

        # Set translate options
        to = gdal.TranslateOptions(format="GTiff", outputSRS="EPSG:25832", noData=-99.0, outputType=gdal.GDT_Float32, callback=progress_callback, callback_data='.')

        # Add tif name and file extension to directory file path.
        out_tif = os.path.join(out_pth, c.REPROJECTED_TIF25832)

        # Convert vrt file to tif
        gdal.Translate(out_tif, out_vrt, options=to)

        return out_tif

        # Close raster
        out_vrt = None
        out_tif = None


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
            log('EXTENSION NOT ASC OR TIF: Filename: {}', (str(f_name)))
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

    return in_ds

    dst_ds = None
    in_ds = None

def reproject_3035(in_ras, out_ras):
    """
    Reproject the final output raster to EPSG:3035.
    """

    # Add reprojected tif name and file extension to directory file path.
    out_pth_ext = os.path.join(out_ras, c.REPROJECTED_TIF3035)

    # Write progress to console
    log_console("\nRunning GDAL Warp to reproject final EPSG:25832 raster to final EPSG:3035...")

    # Reproject
    gdal.Warp(destNameOrDestDS=out_pth_ext, srcDSOrSrcDSTab=in_ras,
              options=gdal.WarpOptions(format='GTiff', srcSRS='EPSG:25832', dstSRS='EPSG:3035',
                                       outputType=gdal.GDT_Float32, callback=progress_callback, callback_data='.'))

    return out_pth_ext

    out_pth_ext = None

def delete_temp_directory(date_time):
    """
    Delete the temporary sub-directory containing all intermediate files created.
    """

    # Delete temp sub directory
    temp_dir_pth = c.root_temp_dir + "/noise_" + date_time
    shutil.rmtree(temp_dir_pth)
