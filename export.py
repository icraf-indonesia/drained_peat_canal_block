import os
import re
import csv
import rasterio
import time
import yaml
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
import matplotlib.pyplot as plt

def create_output_folder(scenario_name, days, n_blocks, n_iterations):
    """Creates an output folder for the simulation scenario."""

    # Sanitize the scenario name (remove invalid characters)
    folder_name = re.sub(r'[\\/*?:"<>|]', "_", scenario_name)
    folder_name = f"{folder_name}_{days}_{n_blocks}_{n_iterations}"
    output_path = os.path.join("output", folder_name)

    # Check if the folder already exists
    if os.path.exists(output_path):
        print(f"WARNING: Output folder '{output_path}' already exists. Existing files may be overwritten.")
    else:
        os.makedirs(output_path)

    return output_path



def export_timeseries_to_csv(timestep_data, output_folder, scenario_name, n_blocks, days):
    """Exports timestep-based variables to a CSV file.

    Args:
        timestep_data (list): List of dictionaries containing timestep data.
        output_folder (str): Path to the output folder.
        scenario_name (str): Name of the simulation scenario.
        n_blocks (int): Number of blocks used in the simulation.
        days (int): Number of days simulated.
    """
    filename = f"{scenario_name}_timeseries_nblocks{n_blocks}_days{days}.csv"
    filepath = os.path.join(output_folder, filename)

    with open(filepath, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)  # Use csv.writer instead of csv.DictWriter

        # Write header row
        header = ['d', 'P', 'ET', 'dry_peat_volume', 'wt_track_drained', 'wt_track_notdrained', 'avg_wt']
        writer.writerow(header)

        # Write timestep data rows
        for row in timestep_data:
            writer.writerow(row)  # Write each sublist as a row


def export_raster_to_geotiff(raster_array, output_path, template_raster, metadata, nodata_value=None):
    """Exports a raster array to a GeoTIFF file.

    Args:
        raster_array (numpy.ndarray): The raster array to export.
        output_path (str): The full path to the output GeoTIFF file.
        template_raster (str): Path to a template raster for geospatial metadata.
        metadata (dict): A dictionary of metadata to add to the GeoTIFF file.
        nodata_value (optional, float): Specific no data value to set for the raster. Defaults to None.
    """
    with rasterio.open(template_raster) as src:
        profile = src.profile.copy()

    profile.update(
        dtype=raster_array.dtype,
        count=1,
        height=raster_array.shape[0],
        width=raster_array.shape[1],
    )

    if nodata_value is not None:
        profile.update(nodata=nodata_value)

    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(raster_array, 1)
        dst.update_tags(**metadata)

def export_simulation_summary(output_folder, scenario_name, days, n_iterations, n_blocks,
                                i, water_blocked_canals, avg_wt_over_time, cumulative_Vdp, 
                                total_co2_emissions=None): # Add total_co2_emissions as optional argument
    """Exports the overall simulation summary to a CSV file.

    Args:
        # ... (Your existing docstring arguments)
        total_co2_emissions (float, optional): Total annual CO2 emissions (Mg). 
            Defaults to None (not included in summary if None).
    """
    filename = f"{scenario_name}_summary_{days}_{n_blocks}_{n_iterations}.csv"
    filepath = os.path.join(output_folder, filename)

    # Create the CSV file
    with open(filepath, 'w', newline='') as csvfile:
        fieldnames = ['Time', 'Days', 'Iterations', 'Blocks', 'Current Iteration',
                      'Water Blocked Canals (m)', 'Avg WTD (m)', 'Cumulative Vdp (m^3)']
        
        # Add 'Total CO2 Emissions (Mg)' to fieldnames if provided
        if total_co2_emissions is not None:
            fieldnames.append('Total CO2 Emissions (Mg)') 
        
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()  # Write the header row
        
        data_row = {  # Start building the data row
            'Time': time.ctime(),
            'Days': days,
            'Iterations': n_iterations,
            'Blocks': n_blocks,
            'Current Iteration': i,
            'Water Blocked Canals (m)': f"{water_blocked_canals:.2f}",
            'Avg WTD (m)': f"{avg_wt_over_time:.2f}",
            'Cumulative Vdp (m^3)': f"{cumulative_Vdp:.2f}"
        }
        
        # Add total CO2 emissions to the data row if provided
        if total_co2_emissions is not None:
            data_row['Total CO2 Emissions (Mg)'] = f"{total_co2_emissions:.2f}" 

        writer.writerow(data_row) # Write the data row


def export_metadata(output_folder, scenario_name, config, start_time, end_time):
    """Exports simulation metadata to a YAML file.

    Args:
        output_folder (str): Path to the output folder.
        scenario_name (str): Name of the simulation scenario.
        config (dict): The loaded YAML configuration.
        start_time (float): Time at the start of the simulation.
        end_time (float): Time at the end of the simulation.
    """
    elapsed_time_seconds = int(end_time - start_time)  # Convert to integer to remove decimals

    metadata = {
        'scenario_name': scenario_name,
        'date_time': time.ctime(),
        'elapsed_time': f"{elapsed_time_seconds} seconds",  # Add seconds unit
        'parameters': config['hydrology'],
        'data_paths': config['data'],
        'tracking_locations': config['tracking'],
        # Add other relevant metadata as needed
    }

    filename = f"{scenario_name}_metadata.yaml"
    filepath = os.path.join(output_folder, filename)

    with open(filepath, 'w') as f:
        yaml.dump(metadata, f, indent=4)

def create_canal_block_outputs(cr_filepath, hand_picked_dams, output_folder, template_raster, scenario_name, days, n_blocks, n_iterations, output_type='both'):
    """
    Creates a Boolean raster map and/or a point shapefile of canal block locations based on the specified parameters.

    Args:
        cr_filepath (str): Path to the GeoTIFF file containing the canal raster (`cr`).
        hand_picked_dams (list): List of integer values representing canal segments with blocks.
        output_folder (str): The folder to store the output files.
        template_raster (str): Path to a template raster for geospatial metadata.
        scenario_name (str): Name of the simulation scenario.
        days (int): Number of days simulated.
        n_blocks (int): Number of blocks used in the simulation.
        n_iterations (int): Number of iterations.
        output_type (str): Type of output to generate ('raster', 'shapefile', or 'both').
    """

    # Ensure output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # 1. Create Boolean Raster (if required):
    if output_type in ['raster', 'both']:
        with rasterio.open(cr_filepath) as src:
            cr_data = src.read(1)
            profile = src.profile.copy()

        # Create a Boolean raster where selected segments are marked as 1
        canal_block_map = np.isin(cr_data, hand_picked_dams).astype(np.uint8)

        # Update profile for the new raster
        profile.update(
            dtype=canal_block_map.dtype,
            count=1,
        )

        # Define raster output filename and path
        raster_filename = f"canal_blocks_{scenario_name}_{days}_{n_blocks}_{n_iterations}.tif"
        raster_output_path = os.path.join(output_folder, raster_filename)
        
        # Write the Boolean raster to a new file
        with rasterio.open(raster_output_path, 'w', **profile) as dst:
            dst.write(canal_block_map, 1)

    # 2. Create Point Shapefile (if required):
    if output_type in ['shapefile', 'both']:
        # Get coordinates of canal pixels
        with rasterio.open(cr_filepath) as src:
            cr_data = src.read(1)

        # Identify rows and columns of pixels greater than 0
        canal_coords = np.argwhere(cr_data > 0)
        rows, cols = canal_coords[:, 0], canal_coords[:, 1]

        # Open the template raster to get its transform (for converting pixel to map coordinates)
        with rasterio.open(template_raster) as src:
            transform = src.transform

        # Create a list of Point geometries for the selected canal blocks
        points = [Point(transform * (col, row)) for row, col in zip(rows, cols) if cr_data[row, col] in hand_picked_dams]

        # Create a GeoDataFrame with the point geometries
        gdf = gpd.GeoDataFrame({'geometry': points}, crs=src.crs)

        # Define shapefile output filename and path
        shapefile_filename = f"canal_blocks_{scenario_name}_{days}_{n_blocks}_{n_iterations}.shp"
        shapefile_output_path = os.path.join(output_folder, shapefile_filename)
        
        # Write the GeoDataFrame to a shapefile
        gdf.to_file(shapefile_output_path, driver='ESRI Shapefile')


def plot_water_table_response(days, wt_track_drained, wt_track_notdrained, avg_wt_over_time):
    """
    Plots the water table response to drainage over time.

    Args:
        days (int): Number of days for the simulation.
        wt_track_drained (list): Water table depth (WTD) values close to the drained area over time.
        wt_track_notdrained (list): WTD values away from the drained area over time.
        avg_wt_over_time (list or float): Average WTD over time. Can be a list or a single float value.
    """
    plt.figure()

    # Plot WTD close to drained area
    plt.plot(list(range(0, days)), wt_track_drained, label='close to drained')

    # Plot WTD away from drained area
    plt.plot(list(range(0, days)), wt_track_notdrained, label='away from drained')

    # Plot average WTD over time with conditional handling for float or list
    if isinstance(avg_wt_over_time, (int, float)):
        plt.plot(list(range(0, days)), [avg_wt_over_time] * days, label='average',
                 linestyle='--', color='red', alpha=0.5)
    else:
        plt.plot(list(range(0, days)), avg_wt_over_time, label='average',
                 linestyle='--', color='red', alpha=0.5)

    # Set plot labels and title
    plt.xlabel('time (days)')
    plt.ylabel('WTD (m)')
    plt.title('Water Table Response to Drainage Over Time')
    plt.legend()

    # Display the plot
    plt.show()

# Example usage:
# plot_water_table_response(30, wt_track_drained, wt_track_notdrained, avg_wt_over_time)
