#!/bin/bash -e

# Debugging.
set -x

if [[ $(verdi code list -Y localhost -r | wc -l) -eq 0 ]]; then
  echo "Installing QE codes..."
  python -m aiidalab_qe install-qe
else
  echo "Quantum ESPRESSO codes are already installed."
fi