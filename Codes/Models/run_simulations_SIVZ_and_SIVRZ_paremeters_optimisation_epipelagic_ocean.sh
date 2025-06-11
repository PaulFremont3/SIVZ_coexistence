#!/bin/bash

for u in {0..3};
do
	for c in {1..50};
	do
		sbatch run_SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.sbatch $u SIVZ_intra $1 $c
	done;
done;
