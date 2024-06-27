#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=36G                                                        
#SBATCH --time=3:00:00
#SBATCH --job-name=Concatenation
#SBATCH --output=%x_%j.out.txt                                                          
#SBATCH --mail-user=rbarrant@uvm.edu
#SBATCH --mail-type=ALL

export samplesheet=$1
export outdir=$2
echo ${samplesheet} 
echo ${outdir} 

module load singularity/3.7.1
source ~/.bashrc
mamba activate env_nf
export NXF_SINGULARITY_CACHEDIR=/gpfs1/mbsr_tools/NXF_SINGULARITY_CACHEDIR
nextflow run preprocess.nf --samplesheet $1 --outdir $2
nextflow clean -f
mamba deactivate
