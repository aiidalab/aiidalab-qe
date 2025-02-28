#!/bin/bash

set -x

# Disable and rename the default 'localhost' computer to 'localhost-non-hq'.
# This is necessary because we will create a new 'localhost' computer
# configured with HyperQueue as the job management system.
verdi computer disable localhost aiida@localhost
python -c "from aiida import orm, load_profile;load_profile();orm.load_computer(label='localhost').label = 'localhost-non-hq'"


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
