#!/bin/bash

# Untart home file to restore home directory

# check if /home is empty by checking if /home/.do_not_delete exists
if [ ! -f /home/jovyan/.do_not_delete ]; then
  # check if /opt/home.tar exists
  if [ -f /opt/home.tar ]; then
    # untar /opt/home.tar to /home/jovyan to restore home directory
    echo "Untarring /opt/home.tar to /home/jovyan"
    tar -xf /opt/home.tar -C /home/jovyan
    # create a file `/home/jovyan/.do_not_delete` to indicate that home directory has been restored
    touch /home/jovyan/.do_not_delete
  fi
fi

# load the ssh-agent and add the default key generated
# the return code can be non-zero if the ssh-agent is not running
# which will cause the notebook to fail to start so we need to ignore the return code
source /opt/bin/load-singlesshagent.sh || true
