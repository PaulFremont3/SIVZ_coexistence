#!/usr/bin/env bash

# remove file first
if [[ ${4} == '1' ]]; then
	rm "$1_$2_res_optimization_$3_POM.txt"
elif [[ ${4} == '0' ]]; then
	rm "$1_$2_res_optimization_$3.txt"
fi

for u in {1..50}; do
    if [[ ${4} == '1' ]]; then
        cat "results_optimization_params/$1_$2_res_optimization_$3_${u}_POM.txt" >> "$1_$2_res_optimization_$3_POM.txt"
    elif [[ ${4} == '0' ]]; then
        cat "results_optimization_params/$1_$2_res_optimization_$3_${u}.txt" >> "$1_$2_res_optimization_$3.txt"
    fi	
done
