#!/bin/bash

set -x

# Disable and rename the default 'localhost' computer to 'localhost-legacy'.
# This is necessary because we will create a new 'localhost' computer
# configured with HyperQueue as the job management system.
verdi computer disable localhost aiida@localhost
verdi computer relabel localhost localhost-legacy


# computer
verdi computer show ${COMPUTER_LABEL} || verdi computer setup        \
  --non-interactive                                               \
  --label "${COMPUTER_LABEL}"                                        \
  --description "local computer with hyperqueue scheduler"        \
  --hostname "localhost"                                          \
  --transport core.local                                          \
  --scheduler hyperqueue                                          \
  --work-dir /home/${NB_USER}/aiida_run/                          \
  --mpirun-command "mpirun -np {num_cpus}"

verdi computer configure core.local "${COMPUTER_LABEL}"              \
  --non-interactive                                               \
  --safe-interval 5.0
