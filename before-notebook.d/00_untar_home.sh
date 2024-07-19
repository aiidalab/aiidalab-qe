#!/bin/bash
#
home="/home/${NB_USER}"

if [[ ! -d ${home} ]]; then
    echo "Directory $home does not exist!"
    exit 1
fi

# Untar home archive file to restore home directory if it is empty
if [[ $(ls -A ${home} | wc -l) = "0" ]];then
  if [ -f /opt/home.tar.gz ]; then
    echo "Extracting /opt/home.tar to /home/${NB_USER}"
    tar -xf /opt/conda/home.tar -C /home/${NB_USER}
  fi
fi
