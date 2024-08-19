#!/bin/bash

set -x

# XXX: need to make daemon start late
verdi daemon stop || echo "stop fail"

# Setup hyperqueue computer if needed
HQ_COMPUTER="local-hq"

computer_list=$(verdi computer list)
if echo ${computer_list} | grep -q ${HQ_COMPUTER}; then
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
      --mpirun-command "mpirun -np {tot_num_mpiprocs}" 

    verdi computer configure core.local "${HQ_COMPUTER}"              \
      --non-interactive                                               \
      --safe-interval 5.0

    # disable the localhost which is set in base image
    # XXX: this cas also be done before hq computer set and set the hq computer as `localhost`
    verdi computer disable localhost aiida@localhost
fi

verdi daemon start || echo "start fail"
