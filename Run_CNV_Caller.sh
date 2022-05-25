#!/bin/bash
#SBATCH --job-name=CNV_Caller
#SBATCH --mail-user=madyson.colton@sahmri.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=/homes/madyson.colton/logs/CNV_220525.log

# Resources allocation request parameters
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=40960      
#SBATCH --time=48:00:00          
#SBATCH --partition sahmri_cancer_hpc

Rscript cnv_caller.R
