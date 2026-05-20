##
# Unified multi-arch Dockerfile for both AMD64 and ARM64 builds.
# Conditionally installs 'bader' in conda (on x86_64) or compiles from source (on arm64),
# and downloads the correct Hyperqueue (hq) binary for each architecture.
##

###############################################################################
# 1) Global ARGs
###############################################################################
ARG FULL_STACK_VER=2026.1030
ARG UV_VER=0.11.15
ARG QE_VER=7.4
ARG QE_DIR=/opt/conda/envs/quantum-espresso-${QE_VER}
ARG HQ_VER=0.19.0

ARG QE_APP_SRC=/tmp/quantum-espresso
ARG COMPUTER_LABEL="localhost"

ARG HQ_URL_AMD64="https://github.com/It4innovations/hyperqueue/releases/download/v${HQ_VER}/hq-v${HQ_VER}-linux-x64.tar.gz"
ARG HQ_URL_ARM64="https://github.com/It4innovations/hyperqueue/releases/download/v${HQ_VER}/hq-v${HQ_VER}-linux-arm64-linux.tar.gz"
ARG MUON_PKG="aiidalab-qe-muon@git+https://github.com/aiidalab/aiidalab-qe-muon@v1.1.1"
ARG AIIDA_HQ_PKG="aiida-hyperqueue~=0.3.0"

###############################################################################
# 2) uv stage (unchanged)
#    https://docs.astral.sh/uv/guides/integration/docker/
###############################################################################
FROM ghcr.io/astral-sh/uv:${UV_VER} AS uv

###############################################################################
# 3) qe_conda_env stage
#    - Creates the quantum-espresso conda environment in /opt/conda/envs/quantum-espresso-<QE_VER>
#    - Install QuantumEspresso and bader packages
###############################################################################
FROM ghcr.io/aiidalab/full-stack:${FULL_STACK_VER} AS qe_conda_env
ARG QE_VER
ARG QE_DIR
# Docker sets TARGETARCH automatically (e.g. "amd64" or "arm64")
ARG TARGETARCH

USER ${NB_USER}
RUN set -ex; \
    echo "Installing QE plus Bader for x86_64..." && \
    mamba create -p ${QE_DIR} --yes qe=${QE_VER} bader && \
    mamba clean --all -f -y

###############################################################################
# 4) base stage to setup common environment variables
###############################################################################
FROM ghcr.io/aiidalab/full-stack:${FULL_STACK_VER} AS base
ARG AIIDA_HQ_PKG
ARG MUON_PKG
ARG COMPUTER_LABEL

ENV COMPUTER_LABEL=${COMPUTER_LABEL}
ENV AIIDA_HQ_PKG=${AIIDA_HQ_PKG}
ENV MUON_PKG=${MUON_PKG}

# https://docs.astral.sh/uv/reference/environment/
ENV UV_CACHE_DIR=/tmp/uv_cache
ENV UV_CONSTRAINT=${PIP_CONSTRAINT}
ENV UV_LINK_MODE=copy
ENV UV_NO_PROGRESS=1
# Make sure UV installs into the conda environment at /opt/conda
ENV UV_SYSTEM_PYTHON=1

###############################################################################
# 5) build_deps stage
#    - Installs Python dependencies using uv for caching
###############################################################################
FROM base AS build_deps
ARG QE_DIR
ARG QE_APP_SRC

WORKDIR ${QE_APP_SRC}

COPY --chown=${NB_UID}:${NB_GID} src/ ${QE_APP_SRC}/src
COPY --chown=${NB_UID}:${NB_GID} setup.cfg pyproject.toml LICENSE README.md ${QE_APP_SRC}

RUN --mount=from=uv,source=/uv,target=/bin/uv \
    --mount=type=cache,sharing=locked,target=${UV_CACHE_DIR},uid=${NB_UID},gid=${NB_GID} \
    uv pip install --strict .

###############################################################################
# 6) home_build stage
#    - Prepares AiiDA profile, sets up hyperqueue, installs QE codes/pseudos,
#      and archives the home folder (home.tar).
###############################################################################
FROM build_deps AS home_build
ARG QE_DIR
# We'll use these to pick the correct HQ binary
ARG TARGETARCH
ARG HQ_URL_AMD64
ARG HQ_URL_ARM64

#
# Download and unpack the correct hq binary for the architecture:
#
RUN set -ex; \
    if [ "${TARGETARCH}" = "arm64" ]; then \
      wget --no-verbose -c -O hq.tar.gz "${HQ_URL_ARM64}"; \
    else \
      wget --no-verbose -c -O hq.tar.gz "${HQ_URL_AMD64}"; \
    fi && \
    tar xf hq.tar.gz -C /opt/conda/

ENV PSEUDO_FOLDER=/tmp/pseudo

RUN mkdir -p ${PSEUDO_FOLDER} && \
    python -m aiidalab_qe download-pseudos --dest ${PSEUDO_FOLDER}

# NOTE: euphonic must be build separately with --no-build-isolation,
# otherwise it fails to build on arm64 due to missing build-time numpy dependency.
RUN --mount=from=uv,source=/uv,target=/bin/uv \
    --mount=type=cache,sharing=locked,target=${UV_CACHE_DIR},uid=${NB_UID},gid=${NB_GID} \
    uv pip install --strict --no-build-isolation euphonic==1.3.2 && \
    uv pip install --strict ${AIIDA_HQ_PKG} ${MUON_PKG} aiidalab-qe-vibroscopy aiida-bader

COPY ./before-notebook.d/* /usr/local/bin/before-notebook.d/

# TODO: Remove PGSQL and daemon log files, and other unneeded files
RUN --mount=from=qe_conda_env,source=${QE_DIR},target=${QE_DIR} \
    bash /usr/local/bin/before-notebook.d/20_start-postgresql.sh && \
    bash /usr/local/bin/before-notebook.d/40_prepare-aiida.sh && \
    bash /usr/local/bin/before-notebook.d/42_setup-hq-computer.sh && \
    python -m aiidalab_qe install-qe --computer ${COMPUTER_LABEL} && \
    python -m aiidalab_qe install-pseudos --source ${PSEUDO_FOLDER} && \
    # steup code: pythonjob, bader, wannier90 code
    verdi code create core.code.installed --label python --computer=localhost --default-calc-job-plugin pythonjob.pythonjob --filepath-executable=/opt/conda/bin/python -n && \
    verdi code create core.code.installed --label bader --computer=localhost --default-calc-job-plugin bader.bader --filepath-executable=${QE_DIR}/bin/bader -n && \
    verdi code create core.code.installed --label wannier90 --computer=localhost --default-calc-job-plugin wannier90.wannier90 --filepath-executable=/opt/conda/bin/wannier90.x -n && \
    # run post_install for plugin
    python -m aiida_bader post-install && \
    python -m aiidalab_qe_vibroscopy setup-phonopy && \
    python -m aiidalab_qe_muon setup-python3 && \
    # wannier90 plugin need SSSP 1.1
    aiida-pseudo install sssp -v 1.1 -x PBE && \
    aiida-pseudo install sssp -v 1.1 -x PBEsol && \
    verdi daemon stop && \
    mamba run -n aiida-core-services pg_ctl stop && \
    touch /home/${NB_USER}/.FLAG_HOME_INITIALIZED && \
    # NOTE: The work folder is empty but if included clashes with the work folder in a Renku
    # session whose permissions cannot be changed.
    # For the same permisssion reason, the .conda folder clashes with aiidalab-launch.
    # It is usually safe (and preferable) to let .conda be recreated on the fly each time,
    # because .conda typically just holds local environment information, caches, or references
    # to available environments.
    cd /home/${NB_USER} && tar -cf /opt/conda/home.tar --exclude .cache --exclude work --exclude .conda .

###############################################################################
# 7) Final stage
#    - Installs python dependencies again, copies QE env, compiles wannier90,
#      and sets up the final environment.
###############################################################################
FROM base
ARG QE_DIR
ARG QE_APP_SRC
ARG TARGETARCH

USER root
# Build wannier90 for all arches
RUN set -ex; \
    apt-get update && apt-get install -y --no-install-recommends \
    gfortran libblas-dev liblapack-dev liblapack3 openmpi-bin libopenmpi-dev; \
    git clone --depth=1 https://github.com/wannier-developers/wannier90.git /tmp/wannier90; \
    cd /tmp/wannier90; \
    cp config/make.inc.gfort make.inc; \
    echo "COMMS=mpi" >> make.inc; \
    echo "MPIF90=mpif90" >> make.inc; \
    make -j"$(nproc)" wannier; \
    cp wannier90.x /opt/conda/bin/wannier90.x; \
    apt-get remove --purge -y gfortran libblas-dev liblapack-dev libopenmpi-dev && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/wannier90 /tmp/bader

USER ${NB_USER}
WORKDIR /tmp

# Install common dependencies such as pymatgen and pandas via conda,
# to avoid building them when installing via pip
# TODO: Remove this once it is part of the full-stack image.
RUN mamba install aiida-core.atomic_tools -y && \
    mamba clean --all -f -y

# Install dependencies in the final image.
# Use uv cache from the previous build step
RUN --mount=from=uv,source=/uv,target=/bin/uv \
    --mount=type=cache,sharing=locked,target=${UV_CACHE_DIR},uid=${NB_UID},gid=${NB_GID} \
    --mount=from=build_deps,source=${QE_APP_SRC},target=${QE_APP_SRC},rw \
    uv pip install --strict --no-build-isolation euphonic==1.3.2 && \
    uv pip install --strict --compile-bytecode \
      ${QE_APP_SRC} ${AIIDA_HQ_PKG} ${MUON_PKG} aiidalab-qe-vibroscopy aiida-bader

# copy hq binary
COPY --from=home_build /opt/conda/hq /usr/local/bin/
# Copy the QE conda environment
COPY --from=qe_conda_env ${QE_DIR} ${QE_DIR}

USER root

# We exclude 42_setup-hq-computer.sh file because the computer is already setup, thus it is not needed in the final image.
COPY ./before-notebook.d/00_untar-home.sh ./before-notebook.d/43_start-hq.sh /usr/local/bin/before-notebook.d/

# Remove the content of /home/<user>, but keep the folder itself
RUN find /home/${NB_USER}/ -mindepth 1 -delete

ENV QE_APP_FOLDER=/opt/conda/quantum-espresso
# Remove all untracked files and directories.
RUN git clean -dffx || true
COPY --chown=${NB_UID}:${NB_GID} . ${QE_APP_FOLDER}

ENV HOME_TAR="/opt/home.tar"
COPY --from=home_build /opt/conda/home.tar "$HOME_TAR"

USER ${NB_USER}
WORKDIR "/home/${NB_USER}"
