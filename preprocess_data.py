# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 09:57:45 2018

@author: L1817
"""

import numpy as np
import rasterio
import scipy.sparse
import pandas as pd


def read_params(fn=r"/data/params.xlsx"):
    df = pd.read_excel(fn)
    return df


def peat_depth_map(peat_depth_type_arr):
    peat_depth_arr = np.ones(shape=peat_depth_type_arr.shape)
    # information from excel
    peat_depth_arr[peat_depth_type_arr==1] = 0 # depth in meters.
    peat_depth_arr[peat_depth_type_arr==2] = 2.
    peat_depth_arr[peat_depth_type_arr==3] = 4.
    peat_depth_arr[peat_depth_type_arr==4] = 0.
    peat_depth_arr[peat_depth_type_arr==5] = 0.
    peat_depth_arr[peat_depth_type_arr==6] = 2.
    peat_depth_arr[peat_depth_type_arr==7] = 4.
    peat_depth_arr[peat_depth_type_arr==8] = 8.
    
    return peat_depth_arr
      
def read_preprocess_rasters(can_rst_fn, dem_rst_fn, peat_type_rst_fn, peat_depth_rst_fn):
    """
    Deals with issues specific to each  set of input data rasters.
    Output
    """
    with rasterio.open(dem_rst_fn) as dem:
        dem = dem.read(1)
    with rasterio.open(can_rst_fn) as can:
        can_arr = can.read(1)
    with rasterio.open(peat_type_rst_fn) as pt:
        peat_type_arr = pt.read(1)
    with rasterio.open(peat_depth_rst_fn) as pd:
        peat_depth_arr = pd.read(1)
        
        
    #Some small changes to get mask of canals: 1 where canals exist, 0 otherwise
    can_arr[can_arr < 0.5] = 0
    can_arr[abs(can_arr) > 0.5] = 1
    can_arr = np.array(can_arr, dtype=int)
    
    # Convert from numpy no data to -9999.0
    dem[dem <-10] = -9999.0
    dem[np.where(np.isnan(dem))] = -9999.0
    dem[dem > 1e20] = -9999.0 # just in case
    
    # control nodata values
    peat_type_arr[peat_type_arr < 0] = -1
    # fill some nodata values to get same size as dem
    peat_type_arr[(np.where(dem>0.1) and np.where(peat_type_arr <0.1))] = 1.
    
    # control nodata values
    peat_depth_arr[peat_depth_arr < 0] = -1
    
    peat_depth_arr = peat_depth_map(peat_depth_arr) # translate number keys to depths
    
    # fill some nodata values to get same size as dem
    peat_depth_arr[(np.where(dem>0.1) and np.where(peat_depth_arr <0.1))] = 1.
    
    
    # Eliminate rows and columns full of noData values.
    # dem = dem[7:-7, 5:-15] #old
    # can_arr = can_arr[7:-7, 5:-15]
    # peat_type_arr = peat_type_arr[7:-7, 5:-15]
    # peat_depth_arr = peat_depth_arr[7:-7, 5:-15]
    
    return can_arr, dem, peat_type_arr, peat_depth_arr
    

# Build the adjacency matrix up
def _prop_to_neighbours(pixel_coords, rasterized_canals, dem, threshold=0.0):
    """Given a pixel where a canals exists, return list of canals that it would propagate to.
    Info taken for allowing propagation: DEM height.
    Threshold gives strictness of the propagation condition: if 0, then only strictly increasing water tables are propagated.
    """
    padded_can_arr = np.pad(rasterized_canals,pad_width=1,mode='constant')
    padded_dem = np.pad(dem,pad_width=1,mode='constant')
    pc0 = pixel_coords[0] + 1; pc1 = pixel_coords[1] + 1 # plus one bc of padding
    prop_to = []

    candidate_coords_list = ([(pc0-1,pc1-1), (pc0-1,pc1), (pc0-1,pc1+1),
                              (pc0,pc1-1),               (pc0,pc1+1),
                              (pc0+1,pc1-1), (pc0+1,pc1), (pc0+1,pc1+1)])
    for cc in candidate_coords_list:
        if padded_can_arr[cc] > 0: # pixel corresponds to a canal
            if padded_dem[cc] - padded_dem[pc0,pc1] > -threshold:
                prop_to.append(int(padded_can_arr[cc]))      
    return prop_to

def gen_can_matrix_and_raster_from_raster(can_rst_fn, dem_rst_fn):
    """ Gets canal RASTER FILE and generates adjacency matrix.
        SLIGHTLY MODIFIED VERSION OF FORESTCARBON´S FUNCTION
    
    Input:
        - can_rst_fn, str: path of the raster file that contains the canals.
        - rst_dem_fn, str: path of the raster file that will be used as a template for burning canal raster.
        
    Output:
        - matrix: canal adjacency matrix NumPy array.
        - out_arr: canal raster in NumPy array.
        - can_to_raster_list: list. Pixels corresponding to each canal.
    """

    can_arr, dem, _, _ = read_preprocess_rasters(can_rst_fn, dem_rst_fn, dem_rst_fn, dem_rst_fn)

    
    # Convert labels of canals to 1,2,3...
    aux =can_arr.flatten()
    can_flat_arr = np.array(aux)
    aux_dem = dem.flatten()
    

    counter = 1
    for i, value in enumerate(aux):
        if value > 0 and aux_dem[i] > 0: # if there is a canal in that pixel and if we have a dem value for that pixel. (Missing dem data points are labelled as -9999)
            can_flat_arr[i] = counter
            counter += 1
    rasterized_canals = can_flat_arr.reshape(can_arr.shape) # contains labels and positions of canals
    n_canals = int(rasterized_canals.max() + 1) 
    
    # compute matrix AND dict
    matrix = scipy.sparse.lil_matrix((n_canals, n_canals)) # lil matrix is good to build it incrementally

    c_to_r_list = [0] * n_canals
    for coords, label in np.ndenumerate(rasterized_canals):
#        print(float(coords[0])/float(dem.shape[0])*100.0, " per cent of the reading completed" ) 
        if rasterized_canals[coords] > 0: # if coords correspond to a canal.
            c_to_r_list[int(label)] = coords # label=0 is not a canal. But do not delete it, otherwise everything would be corrido.
            propagated_to = _prop_to_neighbours(coords, rasterized_canals, dem, threshold=0.0) 
            for i in propagated_to:
                matrix[int(label), i] = 1 # adjacency matrix of the directed graph
    
    
    matrix_csr = matrix.tocsr() # compressed row format is more efficient for later
    matrix_csr.eliminate_zeros() # happens in place. Frees disk usage.
    
    return matrix_csr, rasterized_canals, c_to_r_list

   # Activate this to run this module as main         
#if __name__== '__main__':
#    datafolder = r"C:\Users\L1817\Dropbox\PhD\Computation\Indonesia_WaterTable\Winrock\Canal_Block_Data\GIS_files\Stratification_layers"
#    preprocess_datafolder = r"C:\Users\L1817\Dropbox\PhD\Computation\Indonesia_WaterTable\Winrock\preprocess"
#    dem_rst_fn = r"\dem_clip_cean.tif"
#    can_rst_fn = r"\can_rst_clipped.tif"
#    
#    cm, cr, c_to_r_list = gen_can_matrix_and_raster_from_raster(can_rst_fn=preprocess_datafolder+can_rst_fn,
#                                                                dem_rst_fn=datafolder+dem_rst_fn)
#    # Pickle outcome!
##    import pickle
#    pickle_folder = r"C:\Users\L1817\Winrock"
#    with open(pickle_folder + r'\50x50_DEM_and_canals.pkl', 'w') as f:
#        pickle.dump([cm, cr, c_to_r_list], f)

##########################################################




