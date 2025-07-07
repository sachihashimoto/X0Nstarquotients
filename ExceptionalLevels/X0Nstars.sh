#!/bin/bash

# levels
values=(200 216 224 225)

# run magma
run_magma() {
    local N=$1
    magma N:=${N} X0NstarCCsolver.m
}

# loop over all levels in parallel
for N in "${values[@]}"; do
    run_magma $N &
done

# wait until all processes have been terminated
wait

echo "All levels done."