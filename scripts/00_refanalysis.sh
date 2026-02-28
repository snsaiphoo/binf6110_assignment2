#!/bin/bash

cd ..
mkdir -p refanalysis
cd refanalysis

# Get protein files from data
l329="../data/l329.faa"
s288c="../data/s288C.faa"

# Identity threshold
threshold=80

# Creating a database for the reference strain
apptainer exec ../containers/blast.sif makeblastdb \
  -in s288c \
  -dbtype prot \
  -out S288C_db

# BLASTP performed on both strains to get the mapping and evaluate percent identities
apptainer exec ../containers/blast.sif blastp \
  -query l329 \
  -db S288C_db \
  -out L329_vs_S288C.tsv \
  -evalue 1e-5 \
  -outfmt 6 \
  -num_threads 8

# leep top hit per L-329 protein (highest bitscore = best match)
# filter those hits at ≥${threshold}% identity
# count how many pass the cutoff
# divide by total L-329 proteins
