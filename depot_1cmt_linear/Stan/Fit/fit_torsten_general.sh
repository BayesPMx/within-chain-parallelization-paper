#!/bin/bash

njobs=(1 2 4 8 24)
runnumber=(1 2 3 4 5 6 7 8 9 10)
chain=(1 2 3 4)

for i in "${njobs[@]}";
  do
    for j in "${runnumber[@]}";
      do
        for k in "${chain[@]}";
          do

            mpiexec -n ${i} ./depot_1cmt_prop_torsten_general sample num_warmup=500 num_samples=1000 \
            data file=../../Data/stan_data.json output file=../Fits/Torsten_General/torsten_general_${i}_jobs_run_${j}_${k}.csv \
            random seed=123 init=../../Data/Inits/inits_${j}_${k}.json &

          done
        wait
      done
    wait
  done


# for i in "${njobs[@]}";
#   do
#     for j in "${runnumber[@]}";
#       do
#         for k in "${chain[@]}";
#           do
# 
#             mpiexec -n ${i} ./depot_1cmt_prop_torsten_general sample num_warmup=300 num_samples=200 save_warmup=1 \
#             data file=../../Data/stan_data_short.json output file=../Fits_mpi_test/torsten_general_short_${i}_jobs_run_${j}_${k}.csv \
#             random seed=123 init=../../Data/Inits_short/inits_short_${j}_${k}.json &
# 
#           done
#         wait
#       done
#     wait
#   done


