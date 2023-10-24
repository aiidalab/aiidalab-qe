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
     # The command wil first install the dependencies and then the app and restart the daemon in the end.
     # But since the aiida profile not yet exists, the daemon restart will fail but it is not a problem.
     # Because we only need the dependencies to be installed.
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
# XXX: a fix needed from aiida-pseudo
RUN ${CONDA_DIR}/bin/pip uninstall -y aiida-pseudo && \
     ${CONDA_DIR}/bin/pip install -U git+https://github.com/aiidateam/aiida-pseudo.git@feature/165/download-only-no-profile
ENV PSEUDO_FOLDER ${CONDA_DIR}/pseudo
RUN mkdir -p ${PSEUDO_FOLDER} && \
    python -m aiidalab_qe download-pseudos --dest ${PSEUDO_FOLDER}

COPY before-notebook.d/* /usr/local/bin/before-notebook.d/

WORKDIR "/home/${NB_USER}"
