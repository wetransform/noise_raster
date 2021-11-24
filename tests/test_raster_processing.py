import numpy as np
import unittest

from noise_raster.raster_processing_standalone import sum_sound_level_3D


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