# -*- coding: utf-8 -*-
"""
/***************************************************************************
 NoiseRaster
								 A QGIS plugin
 This plugin processes noise rasters
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
							  -------------------
		begin				 : 2021-09-06
		git sha				 : $Format:%H$
		copyright			 : (C) 2021 by wetransform gmbh
		email				 : kl@wetransform.to
 ***************************************************************************/

/***************************************************************************
 *																		   *
 *	 This program is free software; you can redistribute it and/or modify  *
 *	 it under the terms of the GNU General Public License as published by  *
 *	 the Free Software Foundation; either version 2 of the License, or	   *
 *	 (at your option) any later version.								   *
 *																		   *
 ***************************************************************************/
"""
from qgis.PyQt.QtCore import QSettings, QTranslator, QCoreApplication
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QAction, QMessageBox
import os

from .raster_processing import sum_sound_level_3D, merge_rasters, vectorize, check_projection, validate_source_format, check_extent, create_raster, build_virtual_raster, reproject, source_raster_list, create_zero_array, set_nodata_value, reproject_3035, delete_temp_directory, create_temp_directory, start_logging
import noise_raster.constants as c

# Initialize Qt resources from file resources.py
from .resources import *
# Import the code for the dialog
from .noise_raster_dialog import NoiseRasterDialog
import os.path


class NoiseRaster:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

		:param iface: An interface instance that will be passed to this class
			which provides the hook by which you can manipulate the QGIS
			application at run time.
		:type iface: QgsInterface
		"""
        # Save reference to the QGIS interface
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'NoiseRaster_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)
            QCoreApplication.installTranslator(self.translator)

        # Create the dialog with elements (after translation) and keep reference
        self.dlg = NoiseRasterDialog()

        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&Noise Raster')

        # Check if plugin was started the first time in current QGIS session
        # Must be set in initGui() to survive plugin reloads
        self.first_start = None

        # Create log file
        logger = start_logging()

    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

		We implement this ourselves since we do not inherit QObject.

		:param message: String for translation.
		:type message: str, QString

		:returns: Translated version of message.
		:rtype: QString
		"""
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('NoiseRaster', message)

    def add_action(
            self,
            icon_path,
            text,
            callback,
            enabled_flag=True,
            add_to_menu=True,
            add_to_toolbar=True,
            status_tip=None,
            whats_this=None,
            parent=None):
        """Add a toolbar icon to the toolbar.

		:param icon_path: Path to the icon for this action. Can be a resource
			path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
		:type icon_path: str

		:param text: Text that should be shown in menu items for this action.
		:type text: str

		:param callback: Function to be called when the action is triggered.
		:type callback: function

		:param enabled_flag: A flag indicating if the action should be enabled
			by default. Defaults to True.
		:type enabled_flag: bool

		:param add_to_menu: Flag indicating whether the action should also
			be added to the menu. Defaults to True.
		:type add_to_menu: bool

		:param add_to_toolbar: Flag indicating whether the action should also
			be added to the toolbar. Defaults to True.
		:type add_to_toolbar: bool

		:param status_tip: Optional text to show in a popup when mouse pointer
			hovers over the action.
		:type status_tip: str

		:param parent: Parent widget for the new action. Defaults None.
		:type parent: QWidget

		:param whats_this: Optional text to show in the status bar when the
			mouse pointer hovers over the action.

		:returns: The action that was created. Note that the action is also
			added to self.actions list.
		:rtype: QAction
		"""

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            # Adds plugin icon to Plugins toolbar
            self.iface.addToolBarIcon(action)

        if add_to_menu:
            self.iface.addPluginToMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action

    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        icon_path = ':/plugins/noise_raster/icon.png'
        self.add_action(
            icon_path,
            text=self.tr(u'Process a noise raster'),
            callback=self.run,
            parent=self.iface.mainWindow())

        # will be set False in run()
        self.first_start = True

    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&Noise Raster'),
                action)
            self.iface.removeToolBarIcon(action)

    def run(self):
        """Run method that performs all the real work"""
        # Only create GUI ONCE in callback, so that it will only load when the plugin is started
        if self.first_start == True:
            self.first_start = False

            # Populate combobox dropdown
            self.dlg.comboBox.addItems(['Lden', 'Lnight'])

        # Clear content in input fields from previous run
        lineEdit1 = self.dlg.mQgsFileWidget_1.lineEdit()
        lineEdit2 = self.dlg.mQgsFileWidget_2.lineEdit()
        lineEdit3 = self.dlg.mQgsFileWidget_3.lineEdit()
        lineEdit4 = self.dlg.mQgsFileWidget_out.lineEdit()
        lineEdit1.clear()
        lineEdit2.clear()
        lineEdit3.clear()
        lineEdit4.clear()

        # show the dialog
        self.dlg.show()
        # Run the dialog event loop
        result = self.dlg.exec_()
        # See if OK was pressed
        if result:
            # Create temp sub-directory
            temp_dir, date_time = create_temp_directory()
            # Maximum of 3 source folder paths.
            # Get individual source folder paths.
            noisePths1 = self.dlg.mQgsFileWidget_1.filePath()
            noisePths2 = self.dlg.mQgsFileWidget_2.filePath()
            noisePths3 = self.dlg.mQgsFileWidget_3.filePath()

            # Check for minimum of one file path
            raslist = source_raster_list(noisePths1, noisePths2, noisePths3)

            # Check extent of each input raster
            if len(raslist) == 1:
                check_extent(raslist[0])
            elif len(raslist) == 2:
                check_extent(raslist[0])
                check_extent(raslist[1])
            else:
                check_extent(raslist[0])
                check_extent(raslist[1])
                check_extent(raslist[2])

            # Check file extension of source data. Must be asc or tif.
            if len(raslist) == 1:
                validate_source_format(raslist[0])
            elif len(raslist) == 2:
                validate_source_format(raslist[0])
                validate_source_format(raslist[1])
            else:
                validate_source_format(raslist[0])
                validate_source_format(raslist[1])
                validate_source_format(raslist[2])

            # Check for existence of CRS definition of each GTiff input raster. asc is assumed to be EPSG:25832.
            if len(raslist) == 1:
                check_projection(raslist[0])
            elif len(raslist) == 2:
                check_projection(raslist[0])
                check_projection(raslist[1])
            else:
                check_projection(raslist[0])
                check_projection(raslist[1])
                check_projection(raslist[2])

            # Reproject tifs to EPSG:25832 translate ascs to tifs.
            if len(raslist) == 1:
                reprojectlist = reproject(raslist[0], temp_dir)
                reprojectlist = [reprojectlist]
            elif len(raslist) == 2:
                reprojectlist1 = reproject(raslist[0], temp_dir)
                reprojectlist2 = reproject(raslist[1], temp_dir)
                reprojectlist = [reprojectlist1, reprojectlist2]
            else:
                reprojectlist1 = reproject(raslist[0], temp_dir)
                reprojectlist2 = reproject(raslist[1], temp_dir)
                reprojectlist3 = reproject(raslist[2], temp_dir)
                reprojectlist = [reprojectlist1, reprojectlist2, reprojectlist3]

            # Selected output folder path
            out = self.dlg.mQgsFileWidget_out.filePath()

            # If more than one source path, run addition
            # If only one source path, skip addition and run vectorization
            if len(raslist) > 1:

                # Merge all input noise rasters in each list to create one merged raster, per list
                mergedlist = merge_rasters(reprojectlist, temp_dir)

                # Create merged virtual raster with multiple bands for addition
                mergedVRT = build_virtual_raster(mergedlist, temp_dir)

                # Create masked array for addition which accounts for no data values
                zeroData = create_zero_array(mergedVRT)

                # Call calculation function
                sound_sum = sum_sound_level_3D(zeroData)

                # Write energetically added array to raster in GTiff format
                out_energetic_ras = create_raster(sound_sum, mergedVRT, out)

                # Set no data value to -99.0
                out_final_ras = set_nodata_value(out_energetic_ras)

                # Get reclassification table
                selectedTableIndex = self.dlg.comboBox.currentIndex()

                # Vectorize energetically added raster including all noise sources
                vectorize(out_final_ras, out, selectedTableIndex, temp_dir)

                # Reproject energetically added raster to EPSG:3035
                reproject_3035(out_final_ras, out)

            else:

                # Merge all input rasters for a single noise source
                out_merged_ras = merge_rasters(reprojectlist, temp_dir, out)

                # Get reclassification table
                selectedTableIndex = self.dlg.comboBox.currentIndex()

                # Vectorize raster
                vectorize(out_merged_ras, out, selectedTableIndex, temp_dir)

                # Reproject raster to EPSG:3035
                reproject_3035(out_merged_ras, out)

            pass

            # Invoke modal dialog to enable user to check intermediate files in the temp folder before deletion.
            QMessageBox.information(self.iface.mainWindow(), "Debug",
                                    "Script has completed. You can find the temporary results in " + temp_dir + ". Press OK to delete the temporary result directory.")

            # Delete temporary sub directory containing intermediate files created
            delete_temp_directory(date_time)


            # Load raster layer created by the reproject_3035 function in QGIS
            ras_layer = os.path.join(out, c.REPROJECTED_TIF3035)
            self.iface.addRasterLayer(ras_layer, "out")

            # Load vector layer created by the vectorize function in QGIS
            poly_layer = os.path.join(out, c.REPROJECTED_SHP3035)
            self.iface.addVectorLayer(poly_layer, "out", "ogr")

            # Display success message bar in QGIS
            self.iface.messageBar().pushMessage(
                "Success", "Output file written at " + poly_layer,
                level=3, duration=3)