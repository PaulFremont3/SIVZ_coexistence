#!/bin/bash

for u in {0..3};
do
	echo $u;
	for c in {1..50};
	do
		#u=3;
		echo $u;
		sbatch run_SIVZ_and_SIVRZ_parameters_optimisation_epipelagic_ocean.sbatch $u SIVZ_intra $1 $c
	done;
done;
