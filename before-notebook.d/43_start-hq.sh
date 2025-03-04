#!/bin/bash

set -x

# NOTE: this cgroup folder hierachy is based on cgroupv2
# if the container is open in system which has cgroupv1 the image build procedure will fail.
# Since the image is mostly for demo server where we know the machine and OS I supposed
# it should have cgroupv2 (> Kubernetes v1.25).
# We only build the server for demo server so it does not require user to have new cgroup.
# But for developers, please update your cgroup version to v2.
# See: https://kubernetes.io/docs/concepts/architecture/cgroups/#using-cgroupv2

# Default to cgroupv1 paths
CPU_QUOTA_PATH="/sys/fs/cgroup/cpu/cpu.cfs_quota_us"
CPU_PERIOD_PATH="/sys/fs/cgroup/cpu/cpu.cfs_period_us"
MEMORY_LIMIT_PATH="/sys/fs/cgroup/memory/memory.limit_in_bytes"

# Fallback if cgroupv2 paths exist
if [ -f /sys/fs/cgroup/cpu.max ]; then
  CPU_QUOTA_PATH="/sys/fs/cgroup/cpu.max"
  CPU_PERIOD_PATH="/sys/fs/cgroup/cpu.max"
fi

if [ -f /sys/fs/cgroup/memory.max ]; then
  MEMORY_LIMIT_PATH="/sys/fs/cgroup/memory.max"
fi

# Compute memory limit
if [ -f "$MEMORY_LIMIT_PATH" ]; then
  MEMORY_LIMIT=$(cat "$MEMORY_LIMIT_PATH")
  if [ "$MEMORY_LIMIT" = "max" ] || [ "$MEMORY_LIMIT" -eq -1 ]; then
    MEMORY_LIMIT=4096
  else
    MEMORY_LIMIT=$(echo "scale=0; $MEMORY_LIMIT / (1024 * 1024)" | bc)
  fi
else
  MEMORY_LIMIT=4096
fi
echo "Memory Limit: ${MEMORY_LIMIT} MiB"

# Compute CPU limit
if [ -f "$CPU_QUOTA_PATH" ] && [ -f "$CPU_PERIOD_PATH" ]; then
  CPU_LIMIT=$(cat "$CPU_QUOTA_PATH")
  CPU_PERIOD=$(cat "$CPU_PERIOD_PATH")

  if [ "$CPU_LIMIT" != "max" ] && [ "$CPU_PERIOD" -ne 0 ]; then
    CPU_NUMBER=$(echo "scale=2; $CPU_LIMIT / $CPU_PERIOD" | bc)
    CPU_LIMIT=$(echo "$CPU_NUMBER / 1" | bc)  # Round down to integer
  else
    CPU_LIMIT=$(nproc)
  fi
else
  CPU_LIMIT=$(nproc)
fi

# Temporary fix for https://github.com/aiidalab/aiidalab-qe/issues/1193
# Ensure CPU_LIMIT is at least 1
if [ "$CPU_LIMIT" -le 0 ]; then
  CPU_LIMIT=1
fi

echo "Number of CPUs allocated: $CPU_LIMIT"

# Start HQ server and worker
run-one-constantly hq server start 1>$HOME/.hq-stdout 2>$HOME/.hq-stderr &
run-one-constantly hq worker start --cpus=${CPU_LIMIT} --resource "mem=sum(${MEMORY_LIMIT})" --no-detect-resources &

# Reset the default memory_per_machine and default_mpiprocs_per_machine
# c.set_default_mpiprocs_per_machine = ${CPU_LIMIT}
# c.set_default_memery_per_machine = ${MEMORY_LIMIT}

# Same as original localhost set job poll interval to 2.0 secs
# In addition, set default mpiprocs and memor per machine
# TODO: this will be run every time the container start, we need a lock file to prevent it.
job_poll_interval="2.0"
computer_name=${COMPUTER_LABEL}
python -c "
from aiida import load_profile; from aiida.orm import load_computer;
load_profile();
load_computer('${computer_name}').set_minimum_job_poll_interval(${job_poll_interval})
load_computer('${computer_name}').set_default_mpiprocs_per_machine(${CPU_LIMIT})
load_computer('${computer_name}').set_default_memory_per_machine(${MEMORY_LIMIT})
"
