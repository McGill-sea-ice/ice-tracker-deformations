# Deformations From Sea Ice Tracker Data

## Overview

This project aims at computing arctic sea-ice deformations from icetracker data (Sentinel-1 and RCM). Two distinct methods can be used to compute sea-ice deformations. 

###### Method M00

In the first method, data points with Latitude/Longitude coordinates are processed. A Delaunay triangulation is performed on these data points and the RIOPS grid is used to create local cartesian coordinate systems for each triangular data cell. The triangulated and converted data set is then used to compute sea-ice deformations following *Bouchat et al. (2020)*.

###### Method M01

In the second method, data points with X/Y coordinates are processed. After performing a Delaunay triangulation, we compute sea-ice deformations following *Bouchat et al. (2020)*.

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
python cartopy_feature_download.py physical
~~~


[1]: https://github.com/SciTools/cartopy/pull/1833
## Usage

In order to launch a data processing experience, the main module must be executed. Assuming we are in the project folder, we can execute the main module using the command that follows:

```bash
# Activate the virtual environment
conda activate icetrackdefs
# Launch the code
python src/SeaIceDeformation/main.py
# Deactivate the environment when you are done
conda deactivate
```

The user can configure the experience by modifying the definitions of the parameters in the configuration file `src/SeaIceDeformation/namelist.ini`. In particular, the `output_folder` in the `IO` section should be modified to point to a filesystem location where one has write permissions.  

## Documentation

To generate PDF documentation for this project, start by accessing the `docs` folder (assuming we are already in the project folder):

```bash
cd docs
```

Finally, write the following command:

```bash
make SeaIceDeformation_doc.pdf
```

*SeaIceDeformation_Methods.pdf* will be stored in the current directory.

## Input Data Location

The March and April 2020 Sentinel-1 and RCM ice tracker data stored in *RCM_dats_2020_MarApr.tar.gz* and *S1_dats_2020_MarApr.tar.gz* is located under `/space/hall4/sitestore/eccc/crd/ccrp/mib001/jf_icetracker_data`.


