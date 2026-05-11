#!/bin/bash

SNAKEFOODS=$@

for file in $SNAKEFOODS
	do	
		echo "Running $file:"
		cp -f snakefood/snakefood${file}.json snakefood.json
		# snakemake --cores 1 BQC2_PostannoQC
		snakemake --cores 1 BQC1_PostprocQC
	done
