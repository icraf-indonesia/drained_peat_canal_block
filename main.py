# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 13:13:45 2018

@author: L1817
"""

import os
import re
import sys
import time
import pickle
import yaml
import numpy as np
import matplotlib.pyplot as plt

import preprocess_data
import utilities
import hydro
import hydro_utils
import read
import export

# Configure GDAL Data Path
CONDA_ENV_PATH = sys.prefix
os.environ['GDAL_DATA'] = CONDA_ENV_PATH + '/Library/share/gdal'

plt.close("all")

config_file = "data/original_data/default_data_parameters.yml"


"""
Load configuration from YAML file
"""
with open(config_file, 'r') as file:
    config = yaml.safe_load(file)

# Access parameters from the configuration
scenario_name = config['general']['scenario_name']
days = config['general']['days']
n_blocks = config['general']['n_blocks']
n_iterations = config['general']['n_iterations']
hand_made_dams = config['general']['hand_made_dams']

block_height = config['hydrology']['block_height']
canal_water_level = config['hydrology']['canal_water_level']
diri_bc = config['hydrology']['diri_bc']
hini = config['hydrology']['hini']
et = config['hydrology']['ET'] * np.ones(days)  # Create an array for ET
timestep = config['hydrology']['timestep']
kadjust = config['hydrology']['Kadjust']

dem_rst_fn = config['data']['dem_rst_fn']
can_rst_fn = config['data']['can_rst_fn']
peat_depth_rst_fn = config['data']['peat_depth_rst_fn']
rainfall_fn = config['data']['rainfall_fn']

hand_picked_dams = config['dams']['hand_picked_dams']

# Load WTD tracking areas from YAML
wt_track_drained_area = tuple(config['tracking']['wt_track_drained_area'])
wt_track_notdrained_area = tuple(config['tracking']['wt_track_notdrained_area'])

cum_vdp_nodams = 21088.453521509597  # Value of dry peat volume without blocks, without precipitation for 3 days. Normalization.

"""
Read and preprocess data
"""
# Read rasters, build up canal connectivity adjacency matrix
if 'CNM' and 'cr' and 'c_to_r_list' not in globals(): # call only if needed
    cnm, cr, c_to_r_list = preprocess_data.gen_can_matrix_and_raster_from_raster(
        can_rst_fn=can_rst_fn, dem_rst_fn=dem_rst_fn
    )

_, dem, peat_type_arr, peat_depth_arr = preprocess_data.read_preprocess_rasters(
    can_rst_fn, dem_rst_fn, peat_depth_rst_fn, peat_depth_rst_fn
)

# Read precipitation
p = read.read_precipitation(rainfall_fn) # precipitation read from separate historical data
if len(p) < days:
    raise ValueError("Rainfall data is shorter than the specified simulation days.")
p = p[:days] # Take only the required number of days

# Even if maps say peat depth is less than 2 meters, 
# the impermeable bottom is at most at 2m.
# This can potentially break the hydrological simulation if the WTD would go below 2m.
print(">>>>> WARNING, OVERWRITING PEAT DEPTH")
peat_depth_arr[peat_depth_arr < 2.0] = 2.0

# Catchment mask: delimit the study area
catchment_mask = np.ones(shape=dem.shape, dtype=bool)
catchment_mask[np.where(dem < -10)] = False  # -99999.0 is current value of dem for nodata points.

# 'Peel' the dem. Dirichlet BC will be applied at the peel.
boundary_mask = utilities.peel_raster(dem, catchment_mask)

# After peeling, catchment_mask should only be the fruit:
catchment_mask[boundary_mask] = False

# Soil types, soil physical properties, and soil depth:
peat_type_masked = peat_type_arr * catchment_mask
peat_bottom_elevation = - peat_depth_arr * catchment_mask  # meters with respect to dem surface. Should be negative!

# Load peatmap soil types' physical properties dictionary.
# Kadjust is hydraulic conductivity multiplier for sapric peat
h_to_tra_and_c_dict, _ = hydro_utils.peat_map_interp_functions(Kadjust=kadjust)

# Transmissivity and storage are computed as: T(h) = T(h) - T(peat depth).
#  These quantities are the latter
tra_to_cut = hydro_utils.peat_map_h_to_tra(
    soil_type_mask=peat_type_masked, gwt=peat_bottom_elevation, h_to_tra_and_C_dict=h_to_tra_and_c_dict
)
sto_to_cut = hydro_utils.peat_map_h_to_sto(
    soil_type_mask=peat_type_masked, gwt=peat_bottom_elevation, h_to_tra_and_C_dict=h_to_tra_and_c_dict
)
sto_to_cut = sto_to_cut * catchment_mask.ravel()

# Water level in canals and list of pixels in canal network.
srfcanlist = [dem[coords] for coords in c_to_r_list]
n_canals = len(c_to_r_list)

owtcanlist = [x - canal_water_level for x in srfcanlist]

"""
MonteCarlo
"""
for i in range(0, n_iterations):

    if i == 0:  # random block configurations
        dam_location = np.random.randint(1, n_canals, n_blocks).tolist()  # Generate random kvector. 0 is not a good position in c_to_r_list
    else:
        prohibited_node_list = [i for i, _ in enumerate(owtcanlist[1:]) if owtcanlist[1:][i] < wt_canals[1:][i]]      # [1:] is to take the 0th element out of the loop
        candidate_node_list = np.array([e for e in range(1, n_canals) if e not in prohibited_node_list])  # remove 0 from the range of possible canals
        dam_location = np.random.choice(candidate_node_list, size=n_blocks)

    if hand_made_dams:
        # HAND-MADE RULE OF DAM POSITIONS TO ADD:
        dam_location = hand_picked_dams

    wt_canals = utilities.place_dams(owtcanlist, srfcanlist, block_height, dam_location, cnm)
    """
    #########################################
                    HYDROLOGY
    #########################################
    """
    ny, nx = dem.shape
    dx = 1.0
    dy = 1.0  # metres per pixel  (Actually, pixel size is 100m x 100m, so all units have to be converted afterwards)

    boundary_arr = boundary_mask * (dem - diri_bc)  # constant Dirichlet value in the boundaries

    ele = dem * catchment_mask

    # Get a pickled phi solution (not ele-phi!) computed before without blocks, independently,
    # and use it as initial condition to improve convergence time of the new solution
    retrieve_transient_phi_sol_from_pickled = False
    if retrieve_transient_phi_sol_from_pickled:
        with open(r"pickled/transient_phi_sol.pkl", 'r') as file:
            phi_ini = pickle.load(file)
        print("transient phi solution loaded as initial condition")

    else:
        phi_ini = ele + hini  # initial h (gwl) in the compartment.
        phi_ini = phi_ini * catchment_mask

    wt_canal_arr = np.zeros((ny, nx))  # (nx,ny) array with wt canal height in corresponding nodes
    for canaln, coords in enumerate(c_to_r_list):
        if canaln == 0:
            continue  # because c_to_r_list begins at 1
        wt_canal_arr[coords] = wt_canals[canaln]

    (
        dry_peat_volume,
        cumulative_Vdp,
        avg_wt_over_time,
        wt_track_drained,
        wt_track_notdrained,
        avg_wt,
        timestep_data
    ) = hydro.hydrology(
        'transient', nx, ny, dx, dy, days, ele, phi_ini, catchment_mask, 
        wt_canal_arr, boundary_arr, peat_type_mask=peat_type_masked, 
        httd=h_to_tra_and_c_dict, tra_to_cut=tra_to_cut, sto_to_cut=sto_to_cut, 
        track_WT_drained_area=wt_track_drained_area, track_WT_notdrained_area=wt_track_notdrained_area, # Pass coordinates here
        diri_bc=diri_bc, neumann_bc=None, plotOpt=True, remove_ponding_water=True,
        P=p, ET=et, dt=timestep
    )

    water_blocked_canals = sum(np.subtract(wt_canals[1:], owtcanlist[1:]))

    #cum_Vdp_nodams = 21088.453521509597 # Value of dry peat volume without any blocks, without any precipitation for 3 days. Normalization.
    print(
        'dry_peat_volume(%) = ', dry_peat_volume / cum_vdp_nodams * 100.0, '\n',
        'water_blocked_canals = ', water_blocked_canals
    )

    """
    Final printings
    """
    fname = r'output/results_mc_3_cumulative.txt'
    if n_iterations > 0:  # only if big enough number of simulated days
        with open(fname, 'a') as output_file:
            output_file.write(
                "\n" + str(i) + "    " + str(dry_peat_volume) + "    "
                + str(n_blocks) + "    " + str(n_iterations) + "    " + str(days) + "    "
                + str(time.ctime()) + "    " + str(water_blocked_canals)
            )


# Create the output folder
output_folder = export.create_output_folder(scenario_name, days, n_blocks, n_iterations) 

# Export timeseries data
export.export_timeseries_to_csv(timestep_data, output_folder, scenario_name, n_blocks, days) 

"""
Save WTD data if simulating a year
"""
fname = r'output/wtd_year_' + str(n_blocks) + '.txt'

if days > 300:
    with open(fname, 'a') as output_file:
        output_file.write(
            "\n %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n" 
            + str(time.ctime()) + " nblocks = " + str(n_blocks) + " ET = " + str(et[0])
            + '\n' + 'drained notdrained mean'
        )
        for j in range(len(wt_track_drained)):
            output_file.write(
                "\n" + str(wt_track_drained[j]) + " " + str(wt_track_notdrained[j]) + " " 
                + str(avg_wt_over_time[j])
            )

"""
Water Table Response to Drainage Over Time
"""
plt.figure()
plt.plot(list(range(0, days)), wt_track_drained, label='close to drained')
plt.plot(list(range(0, days)), wt_track_notdrained, label='away from drained')

# Conditional plotting for avg_wt_over_time with dashed red line and 50% opacity
if isinstance(avg_wt_over_time, (int, float)):
    plt.plot(list(range(0, days)), [avg_wt_over_time] * days, label='average',
            linestyle='--', color='red', alpha=0.5)
else:
    plt.plot(list(range(0, days)), avg_wt_over_time, label='average',
            linestyle='--', color='red', alpha=0.5)

plt.xlabel('time(days)')
plt.ylabel('WTD (m)')
plt.title('Water Table Response to Drainage Over Time')
plt.legend()
plt.show()

