#!/bin/bash

njobs=(18)
runnumber=(1)
chain=(1 2 3 4)

for i in "${njobs[@]}";
  do
    for j in "${runnumber[@]}";
      do
        for k in "${chain[@]}";
          do
      
            mpiexec -n 18 ./depot_1cmt_prop_torsten_general sample num_warmup=300 num_samples=200 save_warmup=1 \
            data file=../../Data/stan_data_short.json output file=../Fits/torsten_general_short_${i}_jobs_run_${j}_${k}c.csv \
            random seed=123 init=../../Data/Inits/inits_short_${j}_${k}.json &
        
      done
    done
  done