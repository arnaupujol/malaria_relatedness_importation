
# Probabilities of malaria importation in southern Mozambique through parasite genomics and mobility patterns

This repository contains the code used to generate the results of the analysis
of the manuscript Pujol et al., 'Probabilities of malaria importation in southern Mozambique through parasite genomics and mobility patterns', submitted to Nature Communications.

Software requirements:
----------------------
The software has been tested using Python 3.9.16 and R4.2.3 in a Linux Ubuntu 20.04.4 LTS,
OS type 64-bit and GNOME version 3.36.8.

The R packages required are: 
- dcifer
- tidyverse
- moire

To install them, first install R4.2.3 and RStudio. Then, open RStudio and run the following commands (install the remotes package before if needed, with `install.packages("remotes")`): 

`install.packages("tidyverse")`

`install.packages("dcifer")`

`remotes::install_github("EPPIcenter/moire")`

Alternatively, you can run this other command to install dcifer: 

`remotes::install_github("EPPIcenter/dcifer", build_vignettes = TRUE, force = TRUE)`

All the python packages that are required are:
- numpy
- pandas
- matplotlib
- statsmodels
- scipy
- firthlogist
- geopandas
- contextily
- jupyter lab

And the following local packages:
- [stat_tools version 1.0](https://github.com/arnaupujol/stat_tools)
- [spatial_tools version 1.0](https://github.com/arnaupujol/spatial_tools)
- [genomics version 1.0](https://github.com/arnaupujol/genomics)
- [genmoz_pytools](URL)


Installation instructions:
--------------------------
The installation of the packages was done using conda. For instruction to
install conda in you computed, we refer to the
[conda website](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

For the whole installation procedure, open a Unix terminal to run all the
commands described below.
First, we recommend creating a virtual environment with Python3.8
for this specific analysis:
`conda create -n py39 python=3.9 anaconda`
And then activate the environment:
`conda activate py39`

Once activated, install the required packages:

`conda install python jupyter jupyterlab numpy pandas matplotlib statsmodels scipy geopandas`

`conda install conda-forge::contextily`

`pip install firthlogist`

Then, install the packages from the repositories [stat_tools 1.0](https://github.com/arnaupujol/stat_tools), 
[genmoz_pytools 1.0](https://github.com/arnaupujol/genmoz_pytools), [spatial_tools 1.0](https://github.com/arnaupujol/spatial_tools) and [genomics 1.0](https://github.com/arnaupujol/genomics).

For this, clone the repositories or download and uncompress the zip file of the release, 
and once there, run:

`python setup.py install`

All the packages should be correctly installed. The whole process should not
take more than one hour, unless the full anaconda package is installed instead
of miniconda (see [conda website](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)).

How to use the demos:
----------------------------

To test the performance of the scripts, go to the `demos` directory and lunch
Jupyter lab, for example by running `jupyter lab` in a Unix terminal. A data 
directory would be required containing the databases required for the analyses. 
The original data of Pujol et al. TITLE, JOURNAL can be requested 
by email to arnau.pujol@isglobal.org. 

All the Jupyter notebooks should run correctly. As they are by default, they
show the expected output of each cell run. 
