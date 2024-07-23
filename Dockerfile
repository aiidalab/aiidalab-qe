# syntax=docker/dockerfile:1
ARG FULL_STACK_VER=2024.1021
ARG UV_VER=0.2.27
ARG QE_VER=7.2
ARG QE_DIR=/opt/conda/envs/quantum-espresso-${QE_VER}

FROM ghcr.io/astral-sh/uv:0.2.27 AS uv

# STAGE 1
# Install QE into conda environment in /opt/conda
# This step is independent from the others and can be run in parallel.
FROM ghcr.io/aiidalab/full-stack:${FULL_STACK_VER} AS qe_conda_env
ARG QE_VER
ARG QE_DIR
RUN mamba create -p "${QE_DIR}" --yes qe="${QE_VER}" && \
     mamba clean --all -f -y && \
     fix-permissions "${CONDA_DIR}"


# STAGE 2
FROM ghcr.io/aiidalab/full-stack:${FULL_STACK_VER}
ARG QE_VER
ARG QE_DIR
# Copy whole repo and pre-install the dependencies and app to the tmp folder.
# In the before notebook scripts the app will be re-installed by moving it to the app folder.
ENV PREINSTALL_APP_FOLDER=${CONDA_DIR}/aiidalab-qe
COPY --chown=${NB_UID}:${NB_GID} . ${PREINSTALL_APP_FOLDER}

USER ${NB_USER}

# Using uv to speed up installation, per docs:
# https://github.com/astral-sh/uv/blob/main/docs/guides/docker.md#using-uv-temporarily
# Use the same constraint file as PIP
ENV UV_CONSTRAINT=${PIP_CONSTRAINT}
RUN  --mount=from=uv,source=/uv,target=/bin/uv \
     cd ${PREINSTALL_APP_FOLDER} && \
     # Remove all untracked files and directories. For example the setup lock flag file.
     git clean -fx && \
     uv pip install --system --no-cache . && \
     fix-permissions "${CONDA_DIR}"

# Download the QE pseudopotentials to the folder for afterware installation.
ENV PSEUDO_FOLDER=${CONDA_DIR}/pseudo
RUN mkdir -p ${PSEUDO_FOLDER} && \
    python -m aiidalab_qe download-pseudos --dest ${PSEUDO_FOLDER}

COPY --from=qe_conda_env "${QE_DIR}" "${QE_DIR}"
# TODO: Remove this once we get rid of 70_prepare-qe-executable.sh
ENV QE_VERSION="$QE_VER"

COPY before-notebook.d/* /usr/local/bin/before-notebook.d/

WORKDIR "/home/${NB_USER}"
