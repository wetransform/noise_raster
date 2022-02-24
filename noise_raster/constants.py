# -*- coding: utf-8 -*-
"""
/***************************************************************************
 NoiseRaster constants
 ***************************************************************************/
"""
import tempfile
import os

REPROJECTED_TIF3035 = "final_3035.tif"

REPROJECTED_TIF25832 = "final_25832.tif"

REPROJECTED_SHP3035 = "final_3035.shp"

REPROJECTED_SHP25832 = "final_25832.shp"

logger = None

root_temp_dir = tempfile.gettempdir().replace(os.sep, '/')
