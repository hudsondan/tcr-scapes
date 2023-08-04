# TCR-SCAPES: Sequence Clustering Algorithms for Prediction of TCR Epitope Specificity

A python repository to accompany {publication pending}

## Installation

To install this repository, first clone from Github:

```
git clone https://github.com/hudsondan/tcr-scapes.git
```

To install, please create a new conda environment from the .yml file provided:

## Set up environments

```
conda env create -f clustox_conda.yml
conda activate clustox_conda
```

## Install CoNGA package for access to tcrdist_cpp

Assumes that the user has created and activated a conda virtual environment
Please see https://github.com/phbradley/conga for more details

```
pip install scanpy louvain
git clone https://github.com/phbradley/conga.git && cd conga/tcrdist_cpp && make && cd ..
pip install -e .
conda install seaborn scikit-learn statsmodels numba pytables
conda install -c conda-forge python-igraph leidenalg louvain notebook
cd ..
```

This system has been tested on macOS Monterey and linux-arm64

## Run

To run the benchmarking analysis with default parameters from within an acitvated environment
and save the outputs to a timestamped file in the results folder

```
python run.py
```
To run the programme with a single model , for example with a Hamming distance model

```
python run.py -m hamming -s True
```
To amend the input data, chain selection, models, or any other parameters, see the model help documentation
```
python run.py --help
```
Input data should be formatted similarly to data/vdjdb.csv, the minimum requirements are cdr3, v and j gene codes (formatted as cdr.x, v.x, j.x where x= alpha or beta), pairing (alpha, beta or paired) and epitope columns

## Analyse

Plotting and statistical analyses can be found in results_publication/stats/stats_plots.Rmd
This file references the results produced for the accompanying publication, which can be found in results/publication
CDR3 amino acid motif production can be accessed in make_motifs.py. This will require prior installation of MUSCLE
https://www.drive5.com/muscle/manual/install.html
NB: you may need to change the MUSCLE executable within make_motifs.py by uncommenting the relevant line (#30/#31)

## Acknowledgements

Many of the model implementations in this framework have been adapted from https://github.com/svalkiers/clusTCR.
The package also makes of a C++ implementation of tcrdist from https://github.com/phbradley/conga.

## Help

For feedback on model functionality, issues with installation or running the package, please raise an issue on GitHub at https://github.com/hudsondan/ClustOx.git
For anything else, please contact dan.hudson@dtc.ox.ac.uk
