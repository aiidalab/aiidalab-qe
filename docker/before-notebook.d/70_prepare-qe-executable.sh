#!/bin/bash -e

# Debugging.
set -x

# Copy quantum espresso env to user space.
mkdir -p /home/${NB_USER}/.conda/envs
if [ ! -d /home/${NB_USER}/.conda/envs/quantum-espresso-${QE_VERSION} ]; then
  ln -s /opt/conda/envs/quantum-espresso /home/${NB_USER}/.conda/envs/quantum-espresso-${QE_VERSION}
fi
