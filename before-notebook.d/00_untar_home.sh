#!/bin/bash
set -u

home="/home/${NB_USER}"
HOME_TAR="/opt/conda/home.tar"

# Untar home archive file to restore home directory if it is empty
if [[ $(ls -A ${home} | wc -l) = "0" ]]; then
  if [[ ! -f $HOME_TAR ]]; then
    echo "File $HOME_TAR does not exist!"
    exit 1
  fi
  if [[ ! -d ${QE_APP_FOLDER} ]]; then
    echo "Folder $QE_APP_FOLDER does not exist!"
    exit 1
  fi

  echo "Extracting $HOME_TAR to $home"
  tar -xf $HOME_TAR -C $home

  echo "Copying directory '$QE_APP_FOLDER' to '$home/apps/'"
  cp -r "$QE_APP_FOLDER" "$home/apps/"
fi
set +u
