# Drained Peat Canal Block Simulator

Drained Peat Canal Block Simulator is an open-source model that developed based on research by [Urzainki et al (2020)](https://doi.org/10.5194/bg-17-4769-2020). This model simulates water table depth in peatland areas with canal networks to assess the impact of canal blocking on restoration efforts. The simulation results show how blocking canals reduces carbon emissions from peat decomposition.

## Getting Started

### What you know before you start

This script runs on Python 3, while some other models may still use Python 2. Users should have some basic knowledge of Python. If you're new to Python, it's recommended to start with the [official Python tutorial](https://docs.python.org/3/tutorial/). You'll need to know how to read and write data from files.

The script is designed to be easy to use, you can simply change the main parameters to run it. However, if you want to fully understand the script, it's helpful to use a preferred [integrated development environment (IDE)](https://github.com/learn-co-curriculum/your-integrated-development-environment). This way, you can easily run and understand each line of the script.

### Installation

Before the model can be used, there are installation steps that you need to follow. Here's a step by step guide.

1.  Ensure [Python](https://www.python.org/downloads/) is intalled on your system.

2.  Install [miniconda](https://githubminicondacom/SmithsonianWorkshops/CodingInPython/blob/master/Week%200/Installing%20miniconda%20on%20Windows.md) and run the Anaconda prompt.

3.  Create a new environment within conda with some packages needed

```         
conda create -n [name of environment] python=3 fipy rasterio pandas xlrd
```

3.  Activate the new environment created

```         
conda activate [name of environment]
```

4.  To run the model, you can utilize an integrated development environment (IDE) like [Spyder](https://github.com/spyder-ide/spyder). Here's a step-by-step guide to installing Spyder via conda.

```         
conda create -c conda-forge -n spyder-env spyder numpy scipy pandas matplotlib sympy cython
```

### Check when the model is properly installed

To check the installation, launch the IDE (e.g., Spyder) and attempt to run the *main.py* file as provided. The correct installation will resulted the plot below:

![Plot result from model](src/images/plot-after-computation.png)

If an error message occurs, please investigate and refer to the *common error* solution [WILL PREPARED SOON]. If the error persists, attempt to debug the script.

## Data

The script needs the following data to run:

No | Data | Type | Format | Parameter
--- | --- | --- | --- | ---
1 | Elevation map | Raster | .tif | `dem_rst_fn`
2 | Peat canal network map | Raster | .tif | `can_rst_fn`
3 | Peat depth and soil type map | Raster | .tif | `peat_depth_rst_fn`
4 | Daily precipitation of selected period of time | Tabular | .xlsx | `rainfall_fn`
5 | Information about canal block height, initial canal water level, etc | Tabular | .xlsx | `params_fn`

In addition to all the required data, there are several parameters that need to be defined to run the model. These include the duration of the simulation (in days), the total number of simulated canal blocks, and the total number of simulation iterations.

### Data Preparation

There are certain guidelines you need to follow to prepare data for simulation.

1. All raster data should share the same extent, cell size, and column-row number. This is essential for consistent data processing in the model, requiring uniform properties across datasets.

2. Ensure preparation for tabular data aligns with the template references for [`rainfall_fn`](https://github.com/icraf-indonesia/drained_peat_canal_block/blob/main/data/original_data/params.xlsx) and [`params_fn`](https://github.com/icraf-indonesia/drained_peat_canal_block/blob/main/data/original_data/params.xlsx)

## Reference

Urzainki, I., LaurC)n, A., Palviainen, M., Haahti, K., Budiman, A., Basuki, I., Netzer, M. and HC6kkC\$, H., 2020. Canal blocking optimization in restoration of drained peatlands. Biogeosciences, 17(19), pp.4769-4784.
