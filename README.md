# Drained Peat Canal Block Simulator

Drained Peat Canal Block Simulator is an open-source model that developed based on research by [Urzainki et al (2020)](https://doi.org/10.5194/bg-17-4769-2020). This model simulates ground water table depth in peatland areas with canal networks to assess the impact of canal blocking on restoration efforts. The simulation results show how blocking canals reduces carbon emissions from peat decomposition.

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

1.  An elevation raster such as DEM in a raster image
2.  A canal network map as a raster image
3.  A soil type and depth map as a raster image
4.  Daily precipitation data for a year of selected period of time

## Reference

Urzainki, I., LaurC)n, A., Palviainen, M., Haahti, K., Budiman, A., Basuki, I., Netzer, M. and HC6kkC\$, H., 2020. Canal blocking optimization in restoration of drained peatlands. Biogeosciences, 17(19), pp.4769-4784.
