# syntax=docker/dockerfile:1
FROM ghcr.io/astral-sh/uv:0.2.18 AS uv
FROM ghcr.io/aiidalab/full-stack:2024.1021

USER ${NB_USER}

ENV QE_VERSION="7.2"

# 1. Install Quantum Espresso into a conda environment
# (we do this here for better caching during Docker build)
RUN mamba create -p /opt/conda/envs/quantum-espresso --yes \
    qe=${QE_VERSION} && \
    mkdir -p /home/${NB_USER}/.conda/envs && \
    ln -s /opt/conda/envs/quantum-espresso /home/${NB_USER}/.conda/envs/quantum-espresso-${QE_VERSION} && \
    mamba clean --all -f -y

# 2. Copy the whole repo
ENV QE_APP_FOLDER=${AIIDALAB_APPS}/quantum-espresso
COPY --chown=${NB_UID}:${NB_GID} . ${QE_APP_FOLDER}

WORKDIR "${QE_APP_FOLDER}"
# 3.Install python dependencies
# Use uv instead of pip to speed up installation, per docs:
# https://github.com/astral-sh/uv/blob/main/docs/guides/docker.md#using-uv-temporarily
# Use the same constraint file as pip
ENV UV_CONSTRAINT=${PIP_CONSTRAINT}
RUN --mount=from=uv,source=/uv,target=/bin/uv \
    # Remove all untracked files and directories.
    git clean -dffx && \
    uv pip install --system --no-cache .

# 4. Prepare AiiDA profile and localhost computer
# 5. Install the QE pseudopotentials and codes
RUN bash /usr/local/bin/before-notebook.d/20_start-postgresql.sh && \
    bash /usr/local/bin/before-notebook.d/40_prepare-aiida.sh && \
    python -m aiidalab_qe install-qe && \
    python -m aiidalab_qe install-pseudos && \
    verdi daemon stop && \
    mamba run -n aiida-core-services pg_ctl stop

USER root
COPY ./before-notebook.d/* /usr/local/bin/before-notebook.d/
RUN fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

WORKDIR "/home/${NB_USER}"
RUN tar -cf /opt/home.tar .

USER ${NB_USER}
