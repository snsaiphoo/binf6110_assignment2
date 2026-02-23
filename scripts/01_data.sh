#!/bin/bash

# Move to the raw data folder
cd ..
mkdir -p raw_data
cd raw_data

module load apptainer

# Stage 1 - 10 days
apptainer exec ../containers/sra-toolkit.sif prefetch SRR10551665
apptainer exec ../containers/sra-toolkit.sif prefetch SRR10551664
apptainer exec ../containers/sra-toolkit.sif prefetch SRR10551663

# Stage 2 - 45 days
apptainer exec ../containers/sra-toolkit.sif prefetch SRR10551662
apptainer exec ../containers/sra-toolkit.sif prefetch SRR10551661
apptainer exec ../containers/sra-toolkit.sif prefetch SRR10551660

# Stage 3 - 71 days
apptainer exec ../containers/sra-toolkit.sif prefetch SRR10551659
apptainer exec ../containers/sra-toolkit.sif prefetch SRR10551658
apptainer exec ../containers/sra-toolkit.sif prefetch SRR10551657


