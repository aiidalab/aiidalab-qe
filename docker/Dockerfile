# syntax=docker/dockerfile:1
FROM base-image

USER ${NB_USER}

# Copy whole repo and pre-install the dependencies.
ENV APP_FOLDER /home/${NB_USER}/apps/quantum-espresso
COPY --chown=${NB_UID}:${NB_GID} --from=src . ${APP_FOLDER}

RUN cd ${APP_FOLDER} && \
     # Remove all untracked files and directories. For example the setup lock flag file.
     git clean -fx && \
     pip install . --no-cache-dir && \
     fix-permissions "${CONDA_DIR}" && \
     fix-permissions "/home/${NB_USER}"

# The app version is used for installing the app when first time the container is started.
ARG APP_VERSION
ENV APP_VERSION ${APP_VERSION}

ARG QE_VERSION
ENV QE_VERSION ${QE_VERSION}
RUN mamba create -p /opt/conda/envs/quantum-espresso --yes \
        qe=${QE_VERSION} \
     && mamba clean --all -f -y && \
     fix-permissions "${CONDA_DIR}" && \
     fix-permissions "/home/${NB_USER}"

COPY before-notebook.d/* /usr/local/bin/before-notebook.d/

WORKDIR "/home/${NB_USER}"
