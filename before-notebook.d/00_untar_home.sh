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
    echo "Untarring /opt/home.tar to /home/${NB_USER}"
    tar -xf /opt/home.tar.gz -C /home/${NB_USER}
  fi
# THIS IS ONLY FOR TESTING!
else
    echo "Creating home2!"
    mkdir -p /home/${NB_USER}/home2
    tar -xzf /opt/home.tar.gz -C /home/${NB_USER}/home2
fi
