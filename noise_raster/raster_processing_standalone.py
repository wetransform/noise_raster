import numpy as np
import os

def sum_sound_level_3D(sound_levels: np.array):
    """
    INPUT: array of dimension (m,n,l) - stack of sound levels on a 2D map
    OUTPUT: array of dimension (m,n) - final values of cumulative sound levels on a 2D map

    m: y-nodes on the grid
    n: x-nodes on the grid
    l: no. of sound levels per cell's node

    """

    if not isinstance(sound_levels, np.ndarray): # TODO: with the unittest is redundant
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