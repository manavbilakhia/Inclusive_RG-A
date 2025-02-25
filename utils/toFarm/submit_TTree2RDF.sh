#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name=TTree2RDF
#SBATCH --mail-type=ALL
#SBATCH --mail-user=manavb@jlab.org
#SBATCH --output=/farm_out/%u/%x-%A_%a-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%A_%a-%j-%N.err
#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH --mem-per-cpu=8000
#SBATCH --time=12:00:00

source /u/home/manavb/myenv_clas12.sh

# Navigate to the working directory (if needed)
cd /w/hallb-scshelf2102/clas12/manavb/grad/Inclusive_RG-A/analysis

# Run the executable with the job array index (if needed)
srun ./executable

