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
    # untar /opt/home.tar to /home/jovyan to restore home directory
    echo "Untarring /opt/home.tar to /home/jovyan"
    tar -xf /opt/home.tar.gz -C /home/${NB_USER}
  fi
fi
