#!/bin/bash

cd ..

module load apptainer

names=($(ls -d raw_data/SRR*/ | awk -F'/' '{print $2}'))

for i in "${names[@]}"
do
	echo "$i"

	apptainer exec containers/salmon.sif salmon quant \
	-i salmon_index \
	-l A \
	-r raw_data/${i}.fastq.gz \
	--validateMappings \
	-o salmon_output/${i}
done



