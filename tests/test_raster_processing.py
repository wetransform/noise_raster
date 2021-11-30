import numpy as np
import unittest
import os

from noise_raster.raster_processing_standalone import sum_sound_level_3D, source_raster_list, check_extent

path_to_current_file = os.path.realpath(__file__)
cd = os.path.dirname(path_to_current_file)


def create_test_input_3D_array():
    array1 = np.full((2, 3), 35)
    array2 = np.full((2, 3), 40)
    array3 = np.full((2, 3), 45)
    final_array = np.stack([array1, array2, array3], axis=0)

    return (final_array)

def create_test_output_2D_array():
    array_out = np.full((2, 3), 46.51133105)

    return (array_out)


def test_sum_sound_level_3D():
    # Given
    test_input_3D_array = create_test_input_3D_array()
    test_output_2D_array = create_test_output_2D_array()
    input_shape = test_input_3D_array.shape

    # When
    output = sum_sound_level_3D(test_input_3D_array)

    # Then
    assert isinstance(output, np.ndarray)
    assert output.shape == input_shape[1:3]
    np.allclose(output,test_output_2D_array,equal_nan=False) #returns True if two arrays are element-wise equal within a tolerance (default=1e-05)


def test_source_raster_list():
    # Given
    folderpath1 = os.path.join(cd, "resources", "folderpath1")
    folderpath2 = os.path.join(cd, "resources", "folderpath2")
    folderpath3 = os.path.join(cd, "resources", "folderpath3")

    # When
    out = source_raster_list(folderpath1, folderpath2, folderpath3)

    # Then
    assert isinstance(out, list)
    assert len(out) == 3
    assert len(out[0]) == 3
    assert len(out[1]) == 1
    assert len(out[2]) == 1


def test_check_extent():
    # Given
    testfilepath = 'https://storage.googleapis.com/gcp-public-data-landsat/LC08/01/042/034/LC08_L1TP_042034_20170616_20170629_01_T1/LC08_L1TP_042034_20170616_20170629_01_T1_B4.TIF'
    print(testfilepath)
    # When
    out = check_extent([testfilepath])

    # Then
# #  #   assert out

