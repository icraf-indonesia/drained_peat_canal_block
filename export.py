import os
import re
import csv
import rasterio
import time

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
        writer = csv.writer(csvfile) # Use csv.writer instead of csv.DictWriter

        # Write header row
        header = ['d', 'P', 'ET', 'dry_peat_volume', 'wt_track_drained', 'wt_track_notdrained', 'avg_wt']
        writer.writerow(header)

        # Write timestep data rows
        writer.writerows(timestep_data)
        
def export_raster_to_geotiff(raster_array, output_path, template_raster, metadata):
    """Exports a raster array to a GeoTIFF file.

    Args:
        raster_array (numpy.ndarray): The raster array to export.
        output_path (str): The full path to the output GeoTIFF file.
        template_raster (str): Path to a template raster for geospatial metadata.
        metadata (dict): A dictionary of metadata to add to the GeoTIFF file.
    """
    with rasterio.open(template_raster) as src:
        profile = src.profile.copy()

    profile.update(
        dtype=raster_array.dtype,
        count=1,
        height=raster_array.shape[0],
        width=raster_array.shape[1],
    )

    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(raster_array, 1)
        dst.update_tags(**metadata)


def export_simulation_summary(output_folder, scenario_name, days, n_iterations, n_blocks, 
                                i, water_blocked_canals, avg_wt_over_time, cumulative_Vdp):
    """Exports the overall simulation summary to a CSV file.

    Args:
        # ... (Your existing docstring arguments)
    """
    filename = f"{scenario_name}_summary_{days}_{n_blocks}_{n_iterations}.csv"
    filepath = os.path.join(output_folder, filename)

    # Create the CSV file
    with open(filepath, 'w', newline='') as csvfile:
        fieldnames = ['Time', 'Days', 'Iterations', 'Blocks', 'Current Iteration', 
                      'Water Blocked Canals (m)', 'Avg WTD (m)', 'Cumulative Vdp (m^3)']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader() # Write the header row
        writer.writerow({  # Write the data row
            'Time': time.ctime(),
            'Days': days,
            'Iterations': n_iterations,
            'Blocks': n_blocks,
            'Current Iteration': i,
            'Water Blocked Canals (m)': f"{water_blocked_canals:.2f}", 
            'Avg WTD (m)': f"{avg_wt_over_time:.2f}",
            'Cumulative Vdp (m^3)': f"{cumulative_Vdp:.2f}" 
        })
