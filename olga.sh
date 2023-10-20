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

#for no in 0 500 1000 5000 10000 50000 100000
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
	python run.py -no 100000 -r 1 -cs beta -tcds True -tcdc True -me 1000 -ex "OLGA"
done
conda deactivate
conda deactivate
module purge
echo "DONE"

