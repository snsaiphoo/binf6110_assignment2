#!/bin/bash

# This script will contain the commands needed to obtain the .sif files needed to run the pipeline
# Creates the containers folder to hold the .sif files

cd ..
mkdir -p containers
cd containers

module load apptainer

# SRA Toolkit version 3.2.1
singularity pull docker://quay.io/biocontainers/sra-tools:3.2.1--h4304569_1
mv sra-tools_3.2.1--h4304569_1.sif sra-toolkit.sif

# FASTQC version 0.11.0_cv8
singularity pull docker://biocontainers/fastqc:v0.11.9_cv8
mv fastqc_v0.11.9_cv8.sif fastqc.sif

# Salmon Quantification version 1.10.3
singularity pull docker://combinelab/salmon:1.10.3
mv salmon_1.10.3.sif salmon.sif
