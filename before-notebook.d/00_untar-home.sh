#!/bin/bash
set -eux

home="/home/${NB_USER}"

ls -lrta "$home"
# Untar home archive file to restore home directory if it is empty
if [[ $(ls -A ${home} | wc -l) = 1 ]]; then
  if [[ ! -f $HOME_TAR ]]; then
    echo "File $HOME_TAR does not exist!"
    exit 1
  fi
  if [[ ! -d ${QE_APP_FOLDER} ]]; then
    echo "Folder $QE_APP_FOLDER does not exist!"
    exit 1
  fi

  echo "Extracting $HOME_TAR to $home"
  tar -xf $HOME_TAR -C "$home"

  echo "Copying directory '$QE_APP_FOLDER' to '$AIIDALAB_APPS'"
  cp -r "$QE_APP_FOLDER" "$AIIDALAB_APPS"
else
  echo "$home folder is not empty!"
  ls -lrta "$home"
fi
set +eux
