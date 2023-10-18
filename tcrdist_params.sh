#!/bin/bash
#SBATCH --time=24:00:00 
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --mail-user=dan.hudson@dtc.ox.ac.uk  # adjust this to match your email address
#SBATCH --mail-type=ALL
#SBATCH --output=/well/fernandes/projects/software/slurm/DH/%j.out

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

for r in 12 24 48 96
do
#	for c in 5 10 50 100 500 1000 5000
#	do
#		python run.py -m tcrdist -cs beta -tcdm kmeans -tcdr $r -tcdh $c -ex "tcrdist optimisation" -ds 1000
#		python run.py -m tcrdist3 -cs beta -tcdm kmeans -tcdr $r -tcdh $c -ex "tcrdist optimisation" -ds 1000
#	done
#        
	for h in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1
        do
                python run.py -m tcrdist -cs beta -tcdm DBSCAN -tcdr $r -tcdh $h -ex "tcrdist optimisation" -ds 1000
                python run.py -m tcrdist3 -cs beta -tcdm DBSCAN -tcdr $r -tcdh $h -ex "tcrdist optimisation" -ds 1000
        done
done
conda deactivate
conda deactivate
module purge
echo "DONE"

