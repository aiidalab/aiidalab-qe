#!/bin/bash

set -x

# NOTE: this cgroup folder hierachy is based on cgroupv2
# if the container is open in system which has cgroupv1 the image build procedure will fail.
# Since the image is mostly for demo server where we know the machine and OS I supposed
# it should have cgroupv2 (> Kubernetes v1.25).
# We only build the server for demo server so it does not require user to have new cgroup.
# But for developers, please update your cgroup version to v2.
# See: https://kubernetes.io/docs/concepts/architecture/cgroups/#using-cgroupv2

# computer memory from runtime
MEMORY_LIMIT=$(cat /sys/fs/cgroup/memory.max)

if [ "$MEMORY_LIMIT" = "max" ]; then
  MEMORY_LIMIT=4096
  echo "No memory limit set, use 4GiB"
else
  MEMORY_LIMIT=$(echo "scale=0; $MEMORY_LIMIT / (1024 * 1024)" | bc)
  echo "Memory Limit: ${MEMORY_LIMIT} MiB"
fi

# Compute number of cpus allocated to the container
CPU_LIMIT=$(awk '{print $1}' /sys/fs/cgroup/cpu.max)
CPU_PERIOD=$(awk '{print $2}' /sys/fs/cgroup/cpu.max)

if [ "$CPU_PERIOD" -ne 0 ]; then
  CPU_NUMBER=$(echo "scale=2; $CPU_LIMIT / $CPU_PERIOD" | bc)
  echo "Number of CPUs allocated: $CPU_NUMBER"

  # for HQ setting round to integer number of CPUs, the left are for system tasks
  CPU_LIMIT=$(echo "scale=0; $CPU_LIMIT / $CPU_PERIOD" | bc)
else
  # if no limit (with local OCI without setting cpu limit, use all CPUs)
  CPU_LIMIT=$(nproc)
  echo "No CPU limit set"
fi

# Start hq server with a worker
run-one-constantly hq server start 1>$HOME/.hq-stdout 2>$HOME/.hq-stderr &
run-one-constantly hq worker start --cpus=${CPU_LIMIT} --resource "mem=sum(${MEMORY_LIMIT})" --no-detect-resources &

# Reset the default memory_per_machine and default_mpiprocs_per_machine
# c.set_default_mpiprocs_per_machine = ${CPU_LIMIT}
# c.set_default_memery_per_machine = ${MEMORY_LIMIT}

# Same as original localhost set job poll interval to 2.0 secs
# In addition, set default mpiprocs and memor per machine
# TODO: this will be run every time the container start, we need a lock file to prevent it.
job_poll_interval="2.0"
computer_name=${HQ_COMPUTER}
python -c "
from aiida import load_profile; from aiida.orm import load_computer;
load_profile();
load_computer('${computer_name}').set_minimum_job_poll_interval(${job_poll_interval})
load_computer('${computer_name}').set_default_mpiprocs_per_machine(${CPU_LIMIT})
load_computer('${computer_name}').set_default_memory_per_machine(${MEMORY_LIMIT})
"
