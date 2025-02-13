#!/bin/bash
set -eux

home="/home/${NB_USER}"

# Untar home archive file to restore home directory if it is empty
if [ ! -e $home/.FLAG_HOME_INITIALIZED ]; then
  if [[ ! -f $HOME_TAR ]]; then
    echo "File $HOME_TAR does not exist!"
    exit 1
  fi
  if [[ ! -d ${QE_APP_FOLDER} ]]; then
    echo "Folder $QE_APP_FOLDER does not exist!"
    exit 1
  fi

  echo "Extracting $HOME_TAR to $home"
  # NOTE: the added flags to tar are to address some errors that occur in k8s
  # tar: .: Cannot utime: Operation not permitted -> --touch solves this problem
  # tar: .: Cannot change mode to rwxr-s---: Operation not permitted -> --no-overwrite-dir solves this problem
  # this only refers to the metadata and permissions of existing directories, not their contents
  tar -xf $HOME_TAR -C "$home" --no-overwrite-dir --touch
else
  echo "$home folder is not empty!"
  ls -lrta "$home"
fi

if [ -d $AIIDALAB_APPS/quantum-espresso ]; then
  echo "Quantum ESPRESSO app does exist"
else
  echo "Copying directory '$QE_APP_FOLDER' to '$AIIDALAB_APPS'"
  cp -r "$QE_APP_FOLDER" "$AIIDALAB_APPS"
fi

set +eux
