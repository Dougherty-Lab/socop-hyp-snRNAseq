#!/usr/bin/env bash

#SBATCH --job-name=nf_scalerna
#SBATCH --output=logs/nf_scalerna_%j.out
#SBATCH --error=logs/nf_scalerna_%j.err
#SBATCH --mem=20G
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=s.sarafinovska@wustl.edu


# load system dependencies
eval $(spack load --sh singularityce@3.10.3)
eval $(spack load --sh nextflow@23.04.4)

tmp=$(mktemp -d /tmp/$USER-singularity-XXXXXX)

mkdir local_tmp
mkdir tmp

export NXF_HOME=/scratch/jdlab/s.sarafinovska/.nextflow
export NXF_SINGULARITY_CACHEDIR=singularity_images
export SINGULARITY_TMPDIR=$tmp
export SINGULARITY_CACHEDIR=$tmp


nextflow run ScaleRna-1.5.0-beta2 \
	-profile singularity \
	-params-file runParams.yml \
	-c /ref/jdlab/software/nextflow/wustl_htcf.config \
	--resultDir 240327_MYTsocop \
	--reporting \
	--cellFinder \
	--outDir 240327_MYTsocop-cellFinder-7

