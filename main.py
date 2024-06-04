"""
Simulates peatland hydrology and exports results.
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


start_time = time.time()

config_file = "data/original_data/default_data_parameters_rand_blocks.yml"

# Load Configuration
with open(config_file, 'r') as file:
    config = yaml.safe_load(file)

# General Configuration
scenario_name = config['general']['scenario_name']
days = config['general']['days']
n_blocks = config['general']['n_blocks']
n_iterations = config['general']['n_iterations']
hand_made_dams = config['general']['hand_made_dams']
peat_raster_values = config['general']['peat_raster_values']

# Hydrology Parameters
block_height = config['hydrology']['block_height']
canal_water_level = config['hydrology']['canal_water_level']
diri_bc = config['hydrology']['diri_bc']
hini = config['hydrology']['hini']
et = config['hydrology']['ET'] * np.ones(days)
timestep = config['hydrology']['timestep']
kadjust = config['hydrology']['Kadjust']
dry_period = config['hydrology']['dry_period']

# coefficients for WTD-CO2 emissions relationship
ALPHA = config['co2_em']['ALPHA']
BETA = config['co2_em']['BETA']

# Data File Paths
dem_rst_fn = config['data']['dem_rst_fn']
can_rst_fn = config['data']['can_rst_fn']
peat_depth_rst_fn = config['data']['peat_depth_rst_fn']
rainfall_fn = config['data']['rainfall_fn']

# Dam Placement (Hand-picked)
hand_picked_dams = config['dams']['hand_picked_dams']

# Tracking Locations
wt_track_drained_area = tuple(config['tracking']['wt_track_drained_area'])
wt_track_notdrained_area = tuple(config['tracking']['wt_track_notdrained_area'])
wt_crosssection_y_val = config['tracking']['wt_crosssection_y_val']

# Normalization Value (from manuscript)
cum_vdp_nodams = 21088.453521509597  # Value of dry peat volume without blocks, without precipitation for 3 days.

# Read and Preprocess Data
# Use caching to avoid re-reading if already in memory
if 'CNM' not in globals() or 'cr' not in globals() or 'c_to_r_list' not in globals(): 
    global CNM, cr, c_to_r_list
    CNM, cr, c_to_r_list = preprocess_data.gen_can_matrix_and_raster_from_raster(
        can_rst_fn=can_rst_fn, dem_rst_fn=dem_rst_fn
    )

_, dem, peat_type_arr, peat_depth_arr = preprocess_data.read_preprocess_rasters(
    can_rst_fn, dem_rst_fn, peat_depth_rst_fn, peat_depth_rst_fn
)

p = read.read_precipitation(rainfall_fn)
if len(p) < days:
    raise ValueError("Rainfall data is shorter than the specified simulation days.")
p = p[:days]

print(">>>>> WARNING, OVERWRITING PEAT DEPTH")
peat_depth_arr[peat_depth_arr < 2.0] = 2.0

# Create Masks
catchment_mask = np.ones(shape=dem.shape, dtype=bool)
catchment_mask[np.where(dem < -10)] = False

boundary_mask = utilities.peel_raster(dem, catchment_mask)
catchment_mask[boundary_mask] = False

peat_type_masked = peat_type_arr * catchment_mask
peat_bottom_elevation = - peat_depth_arr * catchment_mask

# --- Calculate peatland area ----
peatland_area_ha = utilities.calculate_peatland_area(peat_type_masked, 100, peat_types=peat_raster_values) 
print(f"Peatland area: {peatland_area_ha:.2f} ha")

# Hydraulic Properties
h_to_tra_and_c_dict, _ = hydro_utils.peat_map_interp_functions(Kadjust=kadjust)

tra_to_cut = hydro_utils.peat_map_h_to_tra(
    soil_type_mask=peat_type_masked, gwt=peat_bottom_elevation, h_to_tra_and_C_dict=h_to_tra_and_c_dict
)
sto_to_cut = hydro_utils.peat_map_h_to_sto(
    soil_type_mask=peat_type_masked, gwt=peat_bottom_elevation, h_to_tra_and_C_dict=h_to_tra_and_c_dict
)
sto_to_cut = sto_to_cut * catchment_mask.ravel()

# Canal Data
srfcanlist = [dem[coords] for coords in c_to_r_list]
n_canals = len(c_to_r_list)
owtcanlist = [x - canal_water_level for x in srfcanlist]

# Create Output Folder
output_folder = export.create_output_folder(scenario_name, days, n_blocks, n_iterations)

# Run Simulation
for i in range(0, n_iterations):

    if i == 0:  # random block configurations
        dam_location = np.random.randint(1, n_canals, n_blocks).tolist()  # Generate random kvector. 0 is not a valid canal in c_to_r_list
    else:
        prohibited_node_list = [j for j, _ in enumerate(owtcanlist[1:]) if owtcanlist[1:][j] < wt_canals[1:][j]]  # [1:] is to take the 0th element out of the loop
        candidate_node_list = np.array([e for e in range(1, n_canals) if e not in prohibited_node_list])  # remove 0 from the range of possible canals
        dam_location = np.random.choice(candidate_node_list, size=n_blocks).tolist()  # Convert to list

    if hand_made_dams:
        # HAND-MADE RULE OF DAM POSITIONS TO ADD:
        dam_location = hand_picked_dams

    wt_canals = utilities.place_dams(owtcanlist, srfcanlist, block_height, dam_location, CNM)

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
        cumulative_dry_peat_volume,
        dry_peat_volume,
        avg_wt_over_time,
        wt_track_drained,
        wt_track_notdrained,
        timestep_data,
        rast_D_before,
        rast_cwl_before,
        rast_dem_before,
        rast_elev_phi_before,
        rast_D_after,
        rast_cwl_after,
        rast_elev_phi_after,
        avg_wtd_nday  # Get n-day average WTD from hydrology function
    ) = hydro.hydrology(
        'transient', nx, ny, dx, dy, days, ele, phi_ini, catchment_mask,
        wt_canal_arr, boundary_arr, peat_type_mask=peat_type_masked,
        httd=h_to_tra_and_c_dict, tra_to_cut=tra_to_cut, sto_to_cut=sto_to_cut,
        track_WT_drained_area=wt_track_drained_area, track_WT_notdrained_area=wt_track_notdrained_area,
        diri_bc=diri_bc, neumann_bc=None, plotOpt=True, remove_ponding_water=True,
        P=p, ET=et, dt=timestep, n_day_avg=dry_period, output_folder=output_folder,
        y_value=wt_crosssection_y_val
    )

    water_blocked_canals = sum(np.subtract(wt_canals[1:], owtcanlist[1:]))

    print(
        'dry_peat_volume(%) = ', dry_peat_volume / cum_vdp_nodams * 100.0, '\n',
        'water_blocked_canals = ', water_blocked_canals
    )

    #  Estimate total CO2 emissions
    avg_wtd_for_emissions = np.mean(avg_wtd_nday)  # Calculate average WTD for emissions (m)
    co2_emissions_per_ha = (- ALPHA * avg_wtd_for_emissions) + BETA # Equation 8, in Mg ha−1 yr−1
    total_co2_emissions = co2_emissions_per_ha * peatland_area_ha # in Mg/yr

    if days < 365:
        print("WARNING: Simulation days less than 365. CO2 emission estimates may be inaccurate.")
    print(f"Estimated CO2 emissions: {total_co2_emissions:.2f} Mg")

    # Plot water table response over time
    export.plot_water_table_response(days, wt_track_drained, wt_track_notdrained, avg_wt_over_time, output_folder)


# Export Outputs
export.export_timeseries_to_csv(timestep_data, output_folder, scenario_name, n_blocks, days)

export.export_simulation_summary(output_folder, scenario_name, days, n_iterations, n_blocks,
                                    i, water_blocked_canals, avg_wt_over_time, cumulative_dry_peat_volume, 
                                    total_co2_emissions, co2_emissions_per_ha)  # Pass the calculated value

end_time = time.time()
export.export_metadata(output_folder, scenario_name, config, start_time, end_time)

cr = cr.astype(np.float32)
filename = f"cr_{scenario_name}_{days}_{n_blocks}_{n_iterations}.tif"
output_path = os.path.join(output_folder, filename)
metadata = {'description': f'Canal raster for scenario {scenario_name}'}
export.export_raster_to_geotiff(cr, output_path, dem_rst_fn, metadata)

# Create canal block outputs (raster and/or shapefile)
if n_blocks > 0 or hand_made_dams:
    export.create_canal_block_outputs(
        cr_filepath=os.path.join(output_folder, f"cr_{scenario_name}_{days}_{n_blocks}_{n_iterations}.tif"),
        hand_picked_dams=dam_location,  # Use the actual dam locations
        output_folder=output_folder,
        template_raster=dem_rst_fn,
        scenario_name=scenario_name,
        days=days,
        n_blocks=n_blocks,
        n_iterations=n_iterations,
        output_type='both'  # Generate both raster and shapefile
    )

# Export raster data
raster_variables = [
    ('rast_D_before', rast_D_before),
    ('rast_cwl_before', rast_cwl_before),
    ('rast_dem_before', rast_dem_before),
    ('rast_elev_phi_before', rast_elev_phi_before),
    ('rast_D_after', rast_D_after),
    ('rast_cwl_after', rast_cwl_after),
    ('rast_elev_phi_after', rast_elev_phi_after),
]
for var_name, raster_array in raster_variables:
    filename = f"{var_name}_{scenario_name}_{days}_{n_blocks}_{n_iterations}.tif"
    output_path = os.path.join(output_folder, filename)
    metadata = {'description': f'{var_name} for scenario {scenario_name}'}
    export.export_raster_to_geotiff(raster_array, output_path, dem_rst_fn, metadata)

