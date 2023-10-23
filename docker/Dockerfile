# syntax=docker/dockerfile:1
FROM base-image

# Copy whole repo and pre-install the dependencies and app to the tmp folder.
# In the before notebook scripts the app will be re-installed by moving it to the app folder.
ENV PREINSTALL_APP_FOLDER ${CONDA_DIR}/aiidalab-qe
COPY --chown=${NB_UID}:${NB_GID} --from=src . ${PREINSTALL_APP_FOLDER}

USER ${NB_USER}

RUN cd ${PREINSTALL_APP_FOLDER} && \
     # Remove all untracked files and directories. For example the setup lock flag file.
     git clean -fx && \
     # It is important to install from `aiidalab install` to mimic the same installation operation from the app store.
     aiidalab install --yes --python ${CONDA_DIR}/bin/python  "quantum-espresso@file://${PREINSTALL_APP_FOLDER}" && \
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

# Download the QE pseudopotentials to the folder for afterware installation.
# XXX: not working because aiida-pseudo install require aiida profile
RUN python -m aiidalab_qe download-pseudos \
     --dest ${CONDA_DIR}

COPY before-notebook.d/* /usr/local/bin/before-notebook.d/

WORKDIR "/home/${NB_USER}"
