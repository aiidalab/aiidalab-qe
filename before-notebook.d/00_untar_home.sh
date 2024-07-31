#!/bin/bash
set -eux

home="/home/${NB_USER}"

# Untar home archive file to restore home directory if it is empty
if [[ $(ls -A ${home} | wc -l) = "0" ]]; then
  if [[ ! -d $TMP_HOME ]]; then
    echo "Folder $TMP_HOME does not exist!"
    exit 1
  fi
  if [[ ! -d ${QE_APP_FOLDER} ]]; then
    echo "Folder $QE_APP_FOLDER does not exist!"
    exit 1
  fi

  echo "Moving $TMP_HOME to $home"
  mv $TMP_HOME/* "$home"

  echo "Copying directory '$QE_APP_FOLDER' to '$AIIDALAB_APPS'"
  cp -r "$QE_APP_FOLDER" "$AIIDALAB_APPS"
else
  echo "$home folder is not empty!"
  ls -lrta "$home"
fi
set +eux
