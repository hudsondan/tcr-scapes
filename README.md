ClustOx: Benchmarking Unsupervised Clustering Models for inference of T Cell Receptor Specificity

A python repository to accompany {publication pending}


Installation:

To install this repository, first clone from Github:
>> git clone www.github.com/hudsondan/ClustOx 

This benchmarking exercise integrates software making use of both Pip and Conda dependencies. To install, please create a new conda environment from the .yml file provided:

# Set up environments

>> conda env create -f clustox.yml
>> conda activate clustox
>> pip install -r requirements.txt

# Install CoNGA package for access to tcrdist_cpp

# Assumes that the user has created and activated a conda virtual environment
# Please see https://github.com/phbradley/conga for more details

!pip3 install scanpy
!git clone https://github.com/phbradley/conga.git && cd conga/tcrdist_cpp && make
!pip install -e .
!conda install seaborn scikit-learn statsmodels numba pytables
!conda install -c conda-forge python-igraph leidenalg louvain notebook
!conda install -c intel tbb # optional
!pip install scanpy
!pip install fastcluster # optional
!conda install pyyaml #optional for using yaml-formatted configuration files for scripts




`