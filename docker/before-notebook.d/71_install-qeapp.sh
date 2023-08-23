#!/bin/bash -e

# Debugging.
set -x

# Install qeapp if it is not already installed.
if aiidalab list | grep -q quantum-espresso; then
    echo "Quantum ESPRESSO app is already installed."
elif [ -z "${APP_VERSION}" ]; then
    echo "This is a development version of the Quantum ESPRESSO app."
else
    echo "Installing Quantum ESPRESSO app."
    aiidalab install --yes quantum-espresso==${APP_VERSION}

    # Remove the repo folder since it is not needed anymore.
    rm -rf /home/${NB_USER}/apps/aiidalab-qe
fi
