#!/bin/bash -e

# Debugging.
set -x

# Copy quantum espresso env to user space.
mkdir -p /home/${NB_USER}/.conda/envs
if [ ! -d /home/${NB_USER}/.conda/envs/quantum-espresso-${QE_VERSION} ]; then
  ln -s /opt/conda/envs/quantum-espresso-${QE_VERSION} /home/${NB_USER}/.conda/envs/quantum-espresso-${QE_VERSION}

  # Install qe so the progress bar not shown in the notebook when first time using app.
  echo "Installing qe."
  python -m aiidalab_qe install-qe
else
  echo "Quantum ESPRESSO app is already installed."
fi
