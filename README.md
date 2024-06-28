# Deformations From Sea Ice Tracker Data

## Overview

This project aims at computing arctic sea-ice deformations from icetracker data (Sentinel-1 and RCM). Two distinct methods can be used to compute sea-ice deformations. 

###### Method M01

In this method, data points with X/Y coordinates are processed. After performing a Delaunay triangulation, we compute sea-ice deformations following *Bouchat et al. (2020)*.

## Installation

Start by cloning the repository:

```bash
# Check if your SSH key is configured on GitHub
ssh -T git@github.com
# Clone the project
git clone git@github.com:McGill-sea-ice/ice-tracker-deformations.git
```

This project uses a [**conda environment**][conda]. Start by accessing the project folder:

[conda]: https://docs.conda.io/en/latest/miniconda.html

```bash
cd ice-tracker-deformations
```

Create and activate the project's environment (this installs the dependencies):

```bash
conda env create -f environment.yaml
conda activate icetrackdefs
```

Install the Cartopy shapefiles (this would be done automatically by Cartopy, but the URL hardcoded in the Cartopy version we used to require is [out of service][1]):
~~~bash
conda activate icetrackdefs
wget -q https://raw.githubusercontent.com/SciTools/cartopy/master/tools/cartopy_feature_download.py -O $CONDA_PREFIX/bin/cartopy_feature_download.py
python $CONDA_PREFIX/bin/cartopy_feature_download.py physical
~~~


[1]: https://github.com/SciTools/cartopy/pull/1833
## Usage

In order to generate a data set, the main module must be executed. Assuming we are in the project folder, we can execute the main module using the following commands:

```bash
# Activate the virtual environment
conda activate icetrackdefs
# Launch the code
python src/SeaIceDeformation/main.py
# Deactivate the environment when you are done
conda deactivate
```

The user can configure the deformation calculations by modifying the definitions of the parameters in the configuration file `src/SeaIceDeformation/namelist.ini`, and the NetCDF analyses can be configured in `src/SatelliteCoverage`. In particular, the `output_folder` in the `IO` section should be modified to point to a filesystem location where one has write permissions. 

The SIDRR.py tool can be used to analyse the data.
(McGill-sea-ice/SIDRRpy


## Testing procedure

in ./data/ run `test.sh`

It will:
 - run `/src/SeaIceDeformation/main.py` with the `namelist.def` settings.
 - run `/data/report_test.py`

A `PASS` should be obtained for all the variables 

## Documentation

We refer to the manuscript submitted to ESSD (ESSD-2024-227).

## Input Data Location

The Data is located on crunch, in /storage/common/


