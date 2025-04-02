#!/bin/bash

set -x

#############################################
# Detect cgroup version
#############################################
IS_CGROUP_V2=0
if [ -f /sys/fs/cgroup/cgroup.controllers ]; then
  IS_CGROUP_V2=1
fi

#############################################
# Generic limit parser
#############################################

# Converts bytes to MiB (returns -1 if input is "max" or too large)
bytes_to_mib() {
  local val=$1
  if [[ "$val" == "max" ]] || [[ "$val" -ge 9223372036854770000 ]]; then
    echo -1
  else
    local mib
    mib=$(awk "BEGIN {v = $val / 1024 / 1024; print (v > 1024000) ? -1 : int(v)}")
    echo "$mib"
  fi
}

# Calculates CPU count from quota and period (rounds down)
calc_cpu_from_quota() {
  local quota=$1
  local period=$2
  if [[ "$quota" == "max" ]] || [[ "$quota" -lt 0 ]] || [[ "$period" -le 0 ]]; then
    echo -1
  else
    awk "BEGIN {printf \"%d\", int($quota / $period)}"
  fi
}

#############################################
# Memory limit
#############################################
read_memory_limit() {
  local path
  if [ "$IS_CGROUP_V2" -eq 1 ]; then
    path="/sys/fs/cgroup/memory.max"
  else
    path="/sys/fs/cgroup/memory/memory.limit_in_bytes"
  fi

  if [ -f "$path" ]; then
    local val
    val=$(<"$path")
    bytes_to_mib "$val"
  else
    echo -1
  fi
}

#############################################
# CPU limit
#############################################
read_cpu_limit() {
  if [ "$IS_CGROUP_V2" -eq 1 ]; then
    local cpu_file="/sys/fs/cgroup/cpu.max"
    if [ -f "$cpu_file" ]; then
      read -r quota period < "$cpu_file"
      calc_cpu_from_quota "$quota" "$period"
    else
      echo -1
    fi
  else
    local quota_file="/sys/fs/cgroup/cpu/cpu.cfs_quota_us"
    local period_file="/sys/fs/cgroup/cpu/cpu.cfs_period_us"
    if [ -f "$quota_file" ] && [ -f "$period_file" ]; then
      local quota period
      quota=$(<"$quota_file")
      period=$(<"$period_file")
      calc_cpu_from_quota "$quota" "$period"
    else
      echo -1
    fi
  fi
}

#############################################
# Main logic
#############################################
MEMORY_LIMIT=$(read_memory_limit)
CPU_LIMIT=$(read_cpu_limit)

# Fallback to system values if needed
if [ "$MEMORY_LIMIT" -le 0 ]; then
  MEMORY_LIMIT=$(python3 -c "import psutil; print(int(psutil.virtual_memory().total / 1024**2))")
fi
if [ "$CPU_LIMIT" -le 0 ]; then
  CPU_LIMIT=$(python3 -c "import os; print(os.cpu_count())")
fi

echo "Detected memory limit:  ${MEMORY_LIMIT} MiB"
echo "Detected CPU limit:     ${CPU_LIMIT}"

#############################################
# Start HQ server and worker
#############################################
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
