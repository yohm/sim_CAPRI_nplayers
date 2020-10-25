#!/bin/bash
#PJM --rsc-list "node=1000"
#PJM --rsc-list "rscunit=rscunit_ft01"
#PJM --rsc-list "rscgrp=eap-large"
#PJM --rsc-list "elapse=05:00:00"
#PJM --mpi "max-proc-per-node=48"
#PJM -s

export OMP_NUM_THREADS=1
mpiexec -stdout-proc ./%/1000R/stdout -stderr-proc ./%/1000R/stderr ../main_evo_fixation_probs_n3.out 10 1.0 0.0001 2.0 96000 50 123456
