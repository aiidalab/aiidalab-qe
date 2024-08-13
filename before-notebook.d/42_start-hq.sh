#!/bin/bash

set -x

# computer memory from runtime
MEMORY_LIMIT=$(cat /sys/fs/cgroup/memory.max)

if [ "$MEMORY_LIMIT" = "max" ]; then
  MEMORY_LIMIT=1024
  echo "No memory limit set"
else
  MEMORY_LIMIT=$(echo "scale=2; $MEMORY_LIMIT / (1024 * 1024)" | bc)
  echo "Memory Limit: ${MEMORY_LIMIT} MB"
fi

# Compute number of cpus allocated to the container
CPU_LIMIT=$(awk '{print $1}' /sys/fs/cgroup/cpu.max)
CPU_PERIOD=$(awk '{print $2}' /sys/fs/cgroup/cpu.max)

if [ "$CPU_PERIOD" -ne 0 ]; then
  CPU_NUMBER=$(echo "scale=2; $CPU_LIMIT / $CPU_PERIOD" | bc)
  echo "Number of CPUs allocated: $CPU_NUMBER"
else
  # if no limit (with local OCI without setting cpu limit, use all CPUs)
  CPU_NUMBER=$(nproc)
  echo "No CPU limit set"
fi

# Start hq server with a worker
run-one-constantly hq server start 1>$HOME/.hq-stdout 2>$HOME/.hq-stderr &
run-one-constantly hq worker start --cpus=${CPU_NUMBER} --resource "mem=sum(${LOCAL_MEM})" --no-detect-resources &

