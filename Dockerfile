# syntax=docker/dockerfile:1
ARG FULL_STACK_VERSION=2024.1021
ARG QE_VERSION=7.2
ARG UV_VERSION=0.2.27

ARG UV_CACHE_DIR=/tmp/uv_cache
ARG QE_DIR=/opt/conda/envs/quantum-espresso-${QE_VERSION}
FROM ghcr.io/astral-sh/uv:${UV_VERSION} AS uv

# STAGE 1
# Install QE into conda environment in /opt/conda
# We need to do this first, otherwise QE gets installed into home folder as part of
# python -m aiidalab_qe install-qe
FROM ghcr.io/aiidalab/full-stack:${FULL_STACK_VERSION} AS qe_conda_env
ARG QE_VERSION
ARG QE_DIR

USER ${NB_USER}
RUN mamba create -p ${QE_DIR} --yes \
    qe=${QE_VERSION} && \
    mamba clean --all -f -y


# STAGE 2
# Install python dependencies needed to run aiidalab_qe CLI commands
FROM ghcr.io/aiidalab/full-stack:${FULL_STACK_VERSION} AS build_deps
ARG UV_CACHE_DIR
ENV QE_APP_FOLDER=/tmp/quantum-espresso
WORKDIR ${QE_APP_FOLDER}

COPY --chown=${NB_UID}:${NB_GID} src/ ${QE_APP_FOLDER}/src
COPY --chown=${NB_UID}:${NB_GID} setup.cfg pyproject.toml *yaml README.md ${QE_APP_FOLDER}

# Use uv instead of pip to speed up installation, per docs:
# https://github.com/astral-sh/uv/blob/main/docs/guides/docker.md#using-uv-temporarily
# Use the same constraint file as pip
ENV UV_CONSTRAINT=${PIP_CONSTRAINT}
RUN --mount=from=uv,source=/uv,target=/bin/uv \
    uv pip install --strict --system --cache-dir=${UV_CACHE_DIR} .


# STAGE 3:
# - Prepare AiiDA profile and localhost computer
# - Install QE codes and pseudopotentials
# - Archive home folder
FROM build_deps AS home_build
ARG QE_DIR

ENV PSEUDO_FOLDER=/tmp/pseudo
RUN mkdir -p ${PSEUDO_FOLDER} && \
    python -m aiidalab_qe download-pseudos --dest ${PSEUDO_FOLDER}

# TODO: Remove PGSQL and daemon log files, and other unneeded files
RUN --mount=from=qe_conda_env,source=${QE_DIR},target=${QE_DIR} \
    bash /usr/local/bin/before-notebook.d/20_start-postgresql.sh && \
    bash /usr/local/bin/before-notebook.d/40_prepare-aiida.sh && \
    python -m aiidalab_qe install-qe && \
    python -m aiidalab_qe install-pseudos --source ${PSEUDO_FOLDER} && \
    verdi daemon stop && \
    mamba run -n aiida-core-services pg_ctl stop && \
    cd /home/${NB_USER} && tar -cf /opt/conda/home.tar .

# STAGE 4 - Final stage
# - Copy QE env environment
# - Copy home folder archive
# - Copy the whole repo content into the container
# - Install python dependencies
# - Remove all content of home folder
FROM ghcr.io/aiidalab/full-stack:${FULL_STACK_VERSION}
ARG QE_DIR
ARG UV_CACHE_DIR
USER ${NB_USER}

COPY --from=qe_conda_env ${QE_DIR} ${QE_DIR}
COPY --from=home_build /opt/conda/home.tar /opt/conda/home.tar

ENV QE_APP_FOLDER=/opt/conda/quantum-espresso
WORKDIR "${QE_APP_FOLDER}"
COPY --chown=${NB_UID}:${NB_GID} . ${QE_APP_FOLDER}
# Remove all untracked files and directories.
RUN git clean -dffx || true

# 3.Install python dependencies
# Use uv cache from the previous build step
ENV UV_CONSTRAINT=${PIP_CONSTRAINT}
RUN --mount=from=uv,source=/uv,target=/bin/uv \
    --mount=from=build_deps,source=${UV_CACHE_DIR},target=${UV_CACHE_DIR},rw \
    uv pip install --strict --system --compile-bytecode --cache-dir=${UV_CACHE_DIR} . && \
    rm -rf build/ src/aiidalab_qe.egg-info/

USER root
COPY ./before-notebook.d/* /usr/local/bin/before-notebook.d/
RUN fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

# REMOVE HOME
RUN find /home/${NB_USER}/ -delete

WORKDIR "/home/${NB_USER}"
USER ${NB_USER}
