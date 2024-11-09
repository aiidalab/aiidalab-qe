# syntax=docker/dockerfile:1
ARG FULL_STACK_VER=2024.1023
ARG UV_VER=0.4.7
ARG QE_VER=7.2
ARG QE_DIR=/opt/conda/envs/quantum-espresso-${QE_VER}
ARG HQ_VER=0.19.0

ARG UV_CACHE_DIR=/tmp/uv_cache
ARG QE_APP_SRC=/tmp/quantum-espresso
ARG HQ_COMPUTER="localhost-hq"

FROM ghcr.io/astral-sh/uv:${UV_VER} AS uv

# STAGE 1
# Install QE into conda environment in /opt/conda
# This step is largely independent from the others and can run in parallel.
# However, it needs to be done before running `python -m aiidalab_qe install-qe`,
# otherwise QE gets installed into ~/.conda folder.
FROM ghcr.io/aiidalab/full-stack:${FULL_STACK_VER} AS qe_conda_env
ARG QE_VER
ARG QE_DIR

USER ${NB_USER}
RUN mamba create -p ${QE_DIR} --yes qe=${QE_VER} && \
    mamba clean --all -f -y

# STAGE 2
# Install python dependencies needed to run aiidalab_qe CLI commands
# uv package cache from this stage is reused in the final stage as well.
FROM ghcr.io/aiidalab/full-stack:${FULL_STACK_VER} AS build_deps
ARG QE_DIR
ARG UV_CACHE_DIR
ARG QE_APP_SRC

WORKDIR ${QE_APP_SRC}
COPY --chown=${NB_UID}:${NB_GID} src/ ${QE_APP_SRC}/src
COPY --chown=${NB_UID}:${NB_GID} setup.cfg pyproject.toml LICENSE README.md ${QE_APP_SRC}

# Use uv instead of pip to speed up installation, per docs:
# https://github.com/astral-sh/uv/blob/main/docs/guides/docker.md#using-uv-temporarily
# Use the same constraint file as pip
ENV UV_CONSTRAINT=${PIP_CONSTRAINT}
RUN --mount=from=uv,source=/uv,target=/bin/uv \
    uv pip install --strict --system --cache-dir=${UV_CACHE_DIR} .

# STAGE 3
# - Prepare AiiDA profile and localhost computer
# - Prepare hq computer using hyperqueue as scheduler
# - Install QE codes and pseudopotentials
# - Archive home folder
FROM build_deps AS home_build
ARG QE_DIR
ARG HQ_VER
ARG HQ_COMPUTER

# Install hq binary
RUN wget -c -O hq.tar.gz https://github.com/It4innovations/hyperqueue/releases/download/v${HQ_VER}/hq-v${HQ_VER}-linux-x64.tar.gz && \
    tar xf hq.tar.gz -C /opt/conda/

ENV PSEUDO_FOLDER=/tmp/pseudo
RUN mkdir -p ${PSEUDO_FOLDER} && \
    python -m aiidalab_qe download-pseudos --dest ${PSEUDO_FOLDER}

ENV UV_CONSTRAINT=${PIP_CONSTRAINT}
# Install the aiida-hyperqueue
# XXX: fix me after release aiida-hyperqueue
RUN --mount=from=uv,source=/uv,target=/bin/uv \
    --mount=from=build_deps,source=${UV_CACHE_DIR},target=${UV_CACHE_DIR},rw \
     uv pip install --system --strict --cache-dir=${UV_CACHE_DIR} \
     "aiida-hyperqueue@git+https://github.com/aiidateam/aiida-hyperqueue"

COPY ./before-notebook.d/* /usr/local/bin/before-notebook.d/

ENV HQ_COMPUTER=$HQ_COMPUTER

# TODO: Remove PGSQL and daemon log files, and other unneeded files
RUN --mount=from=qe_conda_env,source=${QE_DIR},target=${QE_DIR} \
    bash /usr/local/bin/before-notebook.d/20_start-postgresql.sh && \
    bash /usr/local/bin/before-notebook.d/40_prepare-aiida.sh && \
    bash /usr/local/bin/before-notebook.d/42_setup-hq-computer.sh && \
    python -m aiidalab_qe install-qe --computer ${HQ_COMPUTER} && \
    python -m aiidalab_qe install-pseudos --source ${PSEUDO_FOLDER} && \
    verdi daemon stop && \
    mamba run -n aiida-core-services pg_ctl stop && \
    touch /home/${NB_USER}/.FLAG_HOME_INITIALIZED && \
    cd /home/${NB_USER} && tar -cf /opt/conda/home.tar .

# STAGE 3 - Final stage
# - Install python dependencies
# - Copy QE env environment
# - Remove all content of home folder
# - Copy the whole repo content into the container
# - Copy home folder archive
FROM ghcr.io/aiidalab/full-stack:${FULL_STACK_VER}
ARG QE_DIR
ARG QE_APP_SRC
ARG UV_CACHE_DIR
ARG HQ_COMPUTER
USER ${NB_USER}

WORKDIR /tmp
# Install python dependencies
# Use uv cache from the previous build step
# # Install the aiida-hyperqueue
# # XXX: fix me after release aiida-hyperqueue
ENV UV_CONSTRAINT=${PIP_CONSTRAINT}
RUN --mount=from=uv,source=/uv,target=/bin/uv \
    --mount=from=build_deps,source=${UV_CACHE_DIR},target=${UV_CACHE_DIR},rw \
    --mount=from=build_deps,source=${QE_APP_SRC},target=${QE_APP_SRC},rw \
    uv pip install --strict --system --compile-bytecode --cache-dir=${UV_CACHE_DIR} ${QE_APP_SRC} "aiida-hyperqueue@git+https://github.com/aiidateam/aiida-hyperqueue"

# copy hq binary
COPY --from=home_build /opt/conda/hq /usr/local/bin/

COPY --from=qe_conda_env ${QE_DIR} ${QE_DIR}

USER root

COPY ./before-notebook.d/* /usr/local/bin/before-notebook.d/

ENV HQ_COMPUTER=$HQ_COMPUTER

# Remove content of $HOME
# '-mindepth=1' ensures that we do not remove the home directory itself.
RUN find /home/${NB_USER}/ -mindepth 1 -delete

ENV QE_APP_FOLDER=/opt/conda/quantum-espresso
COPY --chown=${NB_UID}:${NB_GID} . ${QE_APP_FOLDER}
# Remove all untracked files and directories.
RUN git clean -dffx || true

ENV HOME_TAR="/opt/home.tar"
COPY --from=home_build /opt/conda/home.tar "$HOME_TAR"

USER ${NB_USER}
WORKDIR "/home/${NB_USER}"
