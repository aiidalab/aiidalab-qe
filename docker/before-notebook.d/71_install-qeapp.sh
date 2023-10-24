#!/bin/bash -e

# Debugging.
set -x

# Install qeapp if it is not already installed.
if aiidalab list | grep -q quantum-espresso; then
    echo "Quantum ESPRESSO app is already installed."
else
    echo "Installing Quantum ESPRESSO app."
    # Install by move the repo folder that is already in the image.
    mv ${PREINSTALL_APP_FOLDER} /home/${NB_USER}/apps/quantum-espresso
fi

# Install the pseudo libraries if not already installed.
# This can be simplified and accelerated once the following PR is merged:
# https://github.com/aiidateam/aiida-pseudo/pull/135
if aiida-pseudo list | grep -q "no pseudo potential families"; then
    echo "Installing pseudo potential families."
    python -m aiidalab_qe install-pseudos --source ${PSEUDO_FOLDER}
else
    echo "Pseudo potential families are already installed."
fi
