#!/bin/bash
set -e
set -u
set -m

QE_VERSION="6.7"
CONDA_ENV_PREFIX="/home/${SYSTEM_USER}/.conda/envs/quantum-espresso-${QE_VERSION}"
echo $CONDA_ENV_PREFIX
exit 0


# Install QE into dedicated conda environment.
if [ ! -d "${CONDA_ENV_PREFIX}" ]
then
  echo -n "Installing QuantumESPRESSO ${QE_VERSION}... "
  conda create --yes --prefix "${CONDA_ENV_PREFIX}" "qe=${QE_VERSION}" 2>&1 >/dev/null && echo "Done."
fi

# Setup code for QE on localhost.
computer_name=localhost

# Setup Quantum ESPRESSO AiiDA codes.
for code_name in pw projwfc dos ;
do
  # Setup code
  verdi code show ${code_name}@${computer_name}>/dev/null || verdi code setup \
      --non-interactive                                             \
      --label ${code_name}                                          \
      --description "${code_name}.x AiiDAlab container."            \
      --input-plugin quantumespresso.${code_name}                   \
      --computer ${computer_name}                                   \
      --remote-abs-path `conda activate quantum-espresso && which ${code_name}.x` \
	  && echo "Setup ${code_name}	OK"
done
