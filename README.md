# Deformations From Sea Ice Tracker Data

## Overview

This project aims at computing arctic sea ice deformations from icetracker data (Sentinel-1 and RCM). 

The RIOPS grid is used to create data quadrilaterals which are in turn used to compute sea ice deformations.

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

## Usage

