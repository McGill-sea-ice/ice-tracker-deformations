# Deformations From Sea Ice Tracker Data

## Overview

This project aims at computing arctic sea ice deformations from icetracker data (Sentinel-1 and RCM). 

A Delaunay triangulation is performed on data points and the RIOPS grid is used to create local cartesian coordinate systems for each triangular data cell. The processed data set is then used to compute sea ice deformations following *Bouchat et al. (2020)*.

## Installation

Start by cloning the repository:

```
ssh -T git@gitlab.science.gc.ca
git clone git@gitlab.science.gc.ca:bdu002/2021_SeaIceDeformation.git
```

This project uses a **virtual environment**. Start by accessing the project folder:

```
cd 2021_SeaIceDeformation
```

Create and activate the project virtual environment:

```
python3 -m venv .venv --without-pip
curl -sS https://bootstrap.pypa.io/get-pip.py -o get-pip.py
source .venv/bin/activate
```

Install the python dependencies on the virtual environment:

```
python -m pip install -r requirements.txt
```




