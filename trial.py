# # -*- coding: utf-8 -*-
# """
# Created on Mon Jan 29 13:46:08 2024

# @author: FMahezs
# """
# import argparse
# import openpyxl
# import pickle
# import numpy as np
# import matplotlib.pyplot as plt
# import time
import os
import preprocess_data


# input data
preprocessed_datafolder = os.path.abspath('./data/Strat4')
dem_rst_fn = os.path.join(preprocessed_datafolder, 'DTM_metres_clip.tif') 
can_rst_fn = os.path.join(preprocessed_datafolder, 'canals_clip.tif')
peat_depth_rst_fn = os.path.join(preprocessed_datafolder, 'Peattypedepth_clip.tif') # peat depth, peat type in the same raster

abs_path_data = os.path.abspath('./data') # Absolute path to data folder needed for Excel file with parameters
params_fn = os.path.join(abs_path_data, 'params.xlsx')


CNM, cr, c_to_r_list = preprocess_data.gen_can_matrix_and_raster_from_raster(can_rst_fn=can_rst_fn, dem_rst_fn=dem_rst_fn)
