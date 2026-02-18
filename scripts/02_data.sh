#!/bin/bash

cd ../raw_data

module load apptainer

apptainer exec ../containers/sra-toolkit.sif fasterq-dump SRR10551665
apptainer exec ../containers/sra-toolkit.sif fasterq-dump SRR10551664
apptainer exec ../containers/sra-toolkit.sif fasterq-dump SRR10551663

apptainer exec ../containers/sra-toolkit.sif fasterq-dump SRR10551662
apptainer exec ../containers/sra-toolkit.sif fasterq-dump SRR10551661
apptainer exec ../containers/sra-toolkit.sif fasterq-dump SRR10551660

apptainer exec ../containers/sra-toolkit.sif fasterq-dump SRR10551659
apptainer exec ../containers/sra-toolkit.sif fasterq-dump SRR10551658
apptainer exec ../containers/sra-toolkit.sif fasterq-dump SRR10551657

gzip *.fastq
