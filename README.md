# nemo_bathymetry
A small package to create a bathymetry file for a NEMO ocean model using an optimal interpolation technique.

## If you want to use this package:
1. Clone this repository: `git clone https://github.com/christophrenkl/nemo_bathymetry.git`
2. Navigate into the directory and creacte Conda environment: `$ conda env create -f environment.yml`
3. Activate Conda environment: `source activate nemo_bathymetry`

## Create a Bathymetry
Follow the steps below to create a bathymetry. Note that the path to the data and values of parameters have to be set in each script.

### Create Background Field
1. `create_background.py`

### Prepare Observations
1. `extract_obs.py`
2. `clean_obs.py`
3. `datum_correction.py`
4. `obs2grid.py`

### Optimal Interpolation of GEBCO and Observations
1. `optimal_interpolation.py`

### Post-Processing
1. `land_sea_mask.py`
2. `min_max_depth.py`
3. `remove_closed_seas.py`


#### Note:
The structure of this project is inspired by [Cookiecutter Data Science](https://drivendata.github.io/cookiecutter-data-science/).
