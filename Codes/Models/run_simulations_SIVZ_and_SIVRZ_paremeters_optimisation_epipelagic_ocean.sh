#!/bin/bash

for u in {0..3};
do
	echo $1
	for c in {1..50};
	do
		echo $c
		#sbatch run_SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.sbatch $u SIVZ $1 $c
		sbatch run_SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.sbatch $u SIVZ_intra $1 $c
		#sbatch run_SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.sbatch $u SIVRZ $1 $c
		#sbatch run_SIVZ_and_SIVRZ_paremeters_optimisation_epipelagic_ocean.sbatch $u SIVRZ_intra $1 $c
	done;
done;
