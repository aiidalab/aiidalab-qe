# syntax=docker/dockerfile:1
FROM aiidalab/qe:amd64-latest

# pg_ctl: cannot be run as root
USER ${NB_USER}
# Run all scripts in the /usr/local/bin/before-notebook.d directory
# to prepare the environment, then stop the services
RUN set -ex && \
    for script in /usr/local/bin/before-notebook.d/*; do \
        /bin/bash "$script"; \
    done && \
    mamba run -n aiida-core-services pg_ctl status -D /home/${NB_USER}/.postgresql && \
    mamba run -n aiida-core-services pg_ctl -w -D /home/${NB_USER}/.postgresql stop && \
    mamba run -n aiida-core-services rabbitmq-server stop && \
    verdi daemon stop

USER root

# cd to /home/${NB_USER}, and tar home folder to /opt/home.tar
# remove all files in /usr/local/bin/before-notebook.d
# Rm all fils in /home/${NB_USER}, including hidden files
RUN cd /home/${NB_USER} && tar -cf /opt/home.tar . && \
    rm -rf /usr/local/bin/before-notebook.d/* && \
    rm -rf /home/${NB_USER}/* /home/${NB_USER}/.[!.]*

# Copy new scripts and make them executable
COPY ./before-notebook.d /usr/local/bin/before-notebook.d
RUN chmod +x /usr/local/bin/before-notebook.d/*

USER ${NB_USER}

WORKDIR "/home/${NB_USER}"
