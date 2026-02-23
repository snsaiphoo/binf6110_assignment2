#!/bin/bash

cd ..

module load apptainer

apptainer exec containers/salmon.sif salmon index \
-t raw_data/rna.fna \
-i salmon_index \
-k 31


