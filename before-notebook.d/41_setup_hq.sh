#!/bin/bash

set -x

# Setup hyperqueue computer if needed
HQ_COMPUTER="local-hq"
LOCALHOST_MPI_PROCS_PER_MACHINE=2

verdi show computer ${HQ_COMPUTER}
if [[ $? -eq 0 ]]; then
  echo "${HQ_COMPUTER} already setup"
else
    # computer
    verdi computer show ${HQ_COMPUTER} || verdi computer setup        \
      --non-interactive                                               \
      --label "${HQ_COMPUTER}"                                        \
      --description "local computer with hyperqueue scheduler"        \
      --hostname "localhost"                                          \
      --transport core.local                                          \
      --scheduler hyperqueue                                          \
      --work-dir /home/${NB_USER}/aiida_run/                          \
      --mpirun-command "mpirun -np {num_cpus}"                        \
      --mpiprocs-per-machine ${LOCALHOST_MPI_PROCS_PER_MACHINE}

    verdi computer configure core.local "${HQ_COMPUTER}"              \
      --non-interactive                                               \
      --safe-interval 5.0
fi
