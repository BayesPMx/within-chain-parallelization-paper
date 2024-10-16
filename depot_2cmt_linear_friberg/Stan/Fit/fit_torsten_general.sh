#!/bin/bash

njobs=(1 4 12 24)
chain=(1 2 3 4)

for i in "${njobs[@]}";
  do
    for k in "${chain[@]}";
      do

        mpiexec -n ${i} ./depot_2cmt_prop_friberg_prop_torsten_general sample num_warmup=500 num_samples=1000 \
        data file=../../Data/stan_data.json output file=../Fits/Torsten_General/torsten_general_${i}_jobs_${k}.csv \
        random seed=1234 init=../../Data/Inits/inits_${k}.json &

      done
    wait
  done

