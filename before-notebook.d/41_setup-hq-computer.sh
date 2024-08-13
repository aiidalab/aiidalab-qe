#!/bin/bash

set -x

# XXX: need to make daemon start late
verdi daemon stop || echo "stop fail"

# Setup hyperqueue computer if needed
HQ_COMPUTER="local-hq"

# XXX: duplicate as 42_??.sh, set in one place as script and reuse?
# Compute number of cpus allocated to the container
CPU_LIMIT=$(awk '{print $1}' /sys/fs/cgroup/cpu.max)
CPU_PERIOD=$(awk '{print $2}' /sys/fs/cgroup/cpu.max)

if [ "$CPU_PERIOD" -ne 0 ]; then
  CPU_NUMBER=$(echo "scale=2; $CPU_LIMIT / $CPU_PERIOD" | bc)
  echo "Number of CPUs allocated: $CPU_NUMBER"

  # for HQ setting round to integer number of CPUs, the left are for system tasks
  HQ_CPU_NUMBER=$(echo "scale=0; $CPU_LIMIT / $CPU_PERIOD" | bc)
else
  # if no limit (with local OCI without setting cpu limit, use all CPUs)
  HQ_CPU_NUMBER=$(nproc)
  echo "No CPU limit set"
fi

computer_list=$(verdi computer list)
if echo ${computer_list} | grep -q ${HQ_COMPUTER}; then
  echo "${HQ_COMPUTER} already setup"
else
    # computer
    # XXX: upbounded mem??
    verdi computer show ${HQ_COMPUTER} || verdi computer setup        \
      --non-interactive                                               \
      --label "${HQ_COMPUTER}"                                        \
      --description "local computer with hyperqueue scheduler"        \
      --hostname "localhost"                                          \
      --transport core.local                                          \
      --scheduler hyperqueue                                          \
      --work-dir /home/${NB_USER}/aiida_run/                          \
      --mpirun-command "mpirun -np {tot_num_mpiprocs}"                        \
      --mpiprocs-per-machine ${HQ_CPU_NUMBER}

    verdi computer configure core.local "${HQ_COMPUTER}"              \
      --non-interactive                                               \
      --safe-interval 5.0
fi

verdi daemon start || echo "start fail"
