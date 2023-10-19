#!/bin/bash
#SBATCH --time=24:00:00 
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --mail-user=dan.hudson@dtc.ox.ac.uk  # adjust this to match your email address
#SBATCH --mail-type=ALL
#SBATCH --output=/well/fernandes/projects/software/slurm/DH/tcr_scapes/%j.out

module purge
module load Python/3.9.6-GCCcore-11.2.0
#module load MUSCLE/5.0.1428-GCCcore-10.3.0 # Loads MUSCLE for sequence alignment 
export PIP_CACHE_DIR="../.pipcache" # Workaround to ensure library does not overload memory
source /well/fernandes/projects/software/miniconda3/bin/activate
#source venv/bin/activate # Activate the local environment
conda activate tcr_scapes

module list
echo "Loaded modules and activated virtual environment"
echo "Running"

for ds in 10 50 100 500 1000
do
	python run.py -r 25 -cs both -tcds True -tcdc True -ds $ds -me $ds -ex "Full benchmarking - tcrdist subset"
#	python run.py -cs both -tcdc True -ds $ds -me $ds -ex "Full benchmarking - all values" 
done
conda deactivate
conda deactivate
module purge
echo "DONE"

