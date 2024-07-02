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

# Install the pseudo:
# If the group exist, the command will skip that library
python -m aiidalab_qe install-pseudos --source ${PSEUDO_FOLDER}
