import os
import re
import csv

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
