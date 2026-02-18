#!/bin/bash

cd ..
mkdir -p fastqc_reports

module load apptainer

apptainer exec containers/fastqc.sif fastqc raw_data/*.fastq.gz -o fastqc_reports/ - t8
