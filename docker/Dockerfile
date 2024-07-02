# syntax=docker/dockerfile:1
FROM ghcr.io/astral-sh/uv:0.2.18 as uv
FROM base-image

# Copy whole repo and pre-install the dependencies and app to the tmp folder.
# In the before notebook scripts the app will be re-installed by moving it to the app folder.
ENV PREINSTALL_APP_FOLDER ${CONDA_DIR}/aiidalab-qe
COPY --chown=${NB_UID}:${NB_GID} --from=src . ${PREINSTALL_APP_FOLDER}

USER ${NB_USER}

# Using uv to speed up installation, per docs:
# https://github.com/astral-sh/uv/blob/main/docs/guides/docker.md#using-uv-temporarily
# Use the same constraint file as PIP
ENV UV_CONSTRAINT ${PIP_CONSTRAINT}
RUN  --mount=from=uv,source=/uv,target=/bin/uv \
     cd ${PREINSTALL_APP_FOLDER} && \
     # Remove all untracked files and directories. For example the setup lock flag file.
     git clean -fx && \
     # It is important to install from `aiidalab install` to mimic the exact installation operation as
     # from the app store.
     # The command wil first install the dependencies from list by parsing setup config files,
     # (for `aiidalab/aiidalab<23.03.2` the `setup.py` should be in the root folder of the app https://github.com/aiidalab/aiidalab/pull/382).
     # and then the app and restart the daemon in the end.
     # But since the aiida profile not yet exists, the daemon restart will fail but it is not a problem.
     # Because we only need the dependencies to be installed.
     # aiidalab install --yes --python ${CONDA_DIR}/bin/python  "quantum-espresso@file://${PREINSTALL_APP_FOLDER}" && \
     # However, have to use `pip install` explicitly because `aiidalab install` call `pip install --user` which will install the app to `/home/${NB_USER}/.local`.
     # It won't cause issue for docker but for k8s deployment the home folder is not bind mounted to the pod and the dependencies won't be found. (see issue in `jupyter/docker-stacks` https://github.com/jupyter/docker-stacks/issues/815)
     uv pip install --system --no-cache . && \
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
ENV PSEUDO_FOLDER ${CONDA_DIR}/pseudo
RUN mkdir -p ${PSEUDO_FOLDER} && \
    python -m aiidalab_qe download-pseudos --dest ${PSEUDO_FOLDER}

COPY before-notebook.d/* /usr/local/bin/before-notebook.d/

WORKDIR "/home/${NB_USER}"
