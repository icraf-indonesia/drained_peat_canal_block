# # -*- coding: utf-8 -*-
# """
# Created on Wed Feb  7 11:20:35 2024

# @author: FMahezs
# """

# ######## Export cr to raster file #########
# import rasterio
# from rasterio.transform import from_origin
# from rasterio.enums import Resampling

# # Define the output raster file path
# output_tif_path = "output/cr_raster.tif"

# # Define pixel size in meters (adjust as needed)
# pixel_size = 100  # for example, assuming a 100m x 100m pixel size

# # Get the shape of the array
# height, width = cr.shape

# # Define the transformation parameters using the bounding box of the array
# transform = from_origin(0, height * pixel_size, pixel_size, pixel_size)

# # Create the output raster file
# with rasterio.open(output_tif_path, 'w', driver='GTiff', height=height, width=width,
#                     count=1, dtype=str(cr.dtype), crs='EPSG:32748', transform=transform) as dst:
#     # Write the data to the raster file
#     dst.write(cr, 1)

# # print(f"Raster exported to: {output_tif_path}")


# ####### Visualize and export the dam location #########
# import matplotlib.pyplot as plt

# def plot_damLoc(wt_canal_arr, damLocation, c_to_r_list, hand_made_dams, hand_picked_dams, N_ITER, DAYS):
#     plt.figure(figsize=(8, 10))
#     plt.imshow(wt_canal_arr, cmap='gray', interpolation='nearest')

#     dam_indices = [c_to_r_list[dam] for dam in damLocation]
#     dam_rows, dam_cols = zip(*[(r, c) for r, c in dam_indices])

#     plt.scatter(dam_cols, dam_rows, marker='.', color='red', label='Dam Points')

#     cb = hand_made_dams
#     nb = len(damLocation)
#     niter = N_ITER
#     ndays = DAYS

#     plt.legend()

#     image_name = f"WTD_CanalBlock{cb}_nblock{nb}_niter{niter}_{ndays}days_dam_dummyX"
#     title_name = f"Dam Location_nblock:{nb}_niter:{niter}_{ndays}days"

#     plt.title(title_name)
#     plt.savefig(image_name, dpi=300)
#     plt.show()

######### Export to excel the results ##########
# import pandas as pd

# df1
# df2

# with pd.ExcelWriter('canal_block_results.xlsx') as writer:
#     df1.to_excel(writer, sheet_name='Sheet_1')
    
    
    
    
    
    
    
    
    
    
    
    
    
    

