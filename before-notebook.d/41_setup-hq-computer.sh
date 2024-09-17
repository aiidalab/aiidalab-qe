#!/bin/bash

set -x

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
      --mpirun-command "mpirun -np {num_cpus}"

    verdi computer configure core.local "${HQ_COMPUTER}"              \
      --non-interactive                                               \
      --safe-interval 5.0

    # disable the localhost which is set in base image
    verdi computer disable localhost aiida@localhost
fi
