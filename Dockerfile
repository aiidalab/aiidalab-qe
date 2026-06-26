##
# Unified multi-arch Dockerfile for both AMD64 and ARM64 builds.
##

###############################################################################
# 1) Global ARGs
###############################################################################
ARG FULL_STACK_VER=2026.1031
ARG QE_VER=7.4
ARG QE_DIR=/opt/conda/envs/quantum-espresso-${QE_VER}
ARG HQ_VER=0.19.0

ARG QE_APP_SRC=/tmp/quantum-espresso
ARG COMPUTER_LABEL="localhost"

ARG HQ_URL_AMD64="https://github.com/It4innovations/hyperqueue/releases/download/v${HQ_VER}/hq-v${HQ_VER}-linux-x64.tar.gz"
ARG HQ_URL_ARM64="https://github.com/It4innovations/hyperqueue/releases/download/v${HQ_VER}/hq-v${HQ_VER}-linux-arm64-linux.tar.gz"
ARG MUON_PKG="aiidalab-qe-muon@git+https://github.com/aiidalab/aiidalab-qe-muon@v1.1.3"
ARG AIIDA_HQ_PKG="aiida-hyperqueue~=0.3.0"


###############################################################################
# 2) qe_conda_env stage
#    - Creates the quantum-espresso conda environment in /opt/conda/envs/quantum-espresso-<QE_VER>
#    - Install QuantumEspresso and bader packages
###############################################################################
FROM ghcr.io/aiidalab/full-stack:${FULL_STACK_VER} AS qe_conda_env
ARG QE_VER
ARG QE_DIR
# Docker sets TARGETARCH automatically (e.g. "amd64" or "arm64")
ARG TARGETARCH

USER ${NB_USER}
RUN echo "Installing QE and Bader..." && \
    mamba create -p ${QE_DIR} --yes qe=${QE_VER} bader && \
    mamba clean --all -f -y

USER root
# Build wannier90 and copy the binary to $QE_DIR conda env
# TODO: Make a conda-forge package for wannier90!
RUN apt-get -q update && \
    apt-get -q install -y --no-install-recommends gfortran libblas-dev liblapack-dev libopenmpi-dev && \
    git clone --depth=1 https://github.com/wannier-developers/wannier90.git /tmp/wannier90 && \
    cd /tmp/wannier90 && \
    cp config/make.inc.gfort make.inc && \
    echo -e "COMMS=mpi\nMPIF90=mpif90" >> make.inc && \
    make -j wannier && \
    cp wannier90.x ${QE_DIR}/bin/wannier90.x

###############################################################################
# 3) base stage to setup common environment variables
#    NOTE: Any ENVs set here will end up in the final image!
###############################################################################
FROM ghcr.io/aiidalab/full-stack:${FULL_STACK_VER} AS base
ARG AIIDA_HQ_PKG
ARG MUON_PKG
ARG COMPUTER_LABEL

ENV COMPUTER_LABEL=${COMPUTER_LABEL}
ENV AIIDA_HQ_PKG=${AIIDA_HQ_PKG}
ENV MUON_PKG=${MUON_PKG}

ENV QE_APP_FOLDER=/opt/conda/quantum-espresso

###############################################################################
# 4) home_build stage
#    - Prepares AiiDA profile, sets up hyperqueue, installs QE codes/pseudos,
#      and archives the home folder (home.tar).
###############################################################################
FROM base AS home_build
ARG QE_APP_SRC
ARG QE_DIR
# We'll use these to pick the correct HQ binary
ARG TARGETARCH
ARG HQ_URL_AMD64
ARG HQ_URL_ARM64

WORKDIR ${QE_APP_SRC}

COPY --chown=${NB_UID}:${NB_GID} src/ ${QE_APP_SRC}/src
COPY --chown=${NB_UID}:${NB_GID} setup.cfg pyproject.toml LICENSE README.md ${QE_APP_SRC}

RUN python -m pip install --no-user --no-cache-dir ${AIIDA_HQ_PKG}

# Download and unpack the correct hq binary for the architecture:
RUN set -ex; \
    if [ "${TARGETARCH}" = "arm64" ]; then \
      wget --no-verbose -c -O hq.tar.gz "${HQ_URL_ARM64}"; \
    else \
      wget --no-verbose -c -O hq.tar.gz "${HQ_URL_AMD64}"; \
    fi && \
    tar xf hq.tar.gz -C /opt/conda/

# Install common dependencies such as pymatgen and pandas via conda,
# to avoid building them when installing via pip
# TODO: Remove this once it is part of the full-stack image.
RUN mamba install aiida-core.atomic_tools -y && \
    mamba clean --all -f -y

# Install the app and its plugins into the user's local Python environment
RUN python -m pip install --user --no-cache-dir . ${MUON_PKG} aiidalab-qe-vibroscopy aiida-bader

ENV PSEUDO_FOLDER=/tmp/pseudo

RUN mkdir -p ${PSEUDO_FOLDER} && \
    python -m aiidalab_qe download-pseudos --dest ${PSEUDO_FOLDER}

COPY ./before-notebook.d/* /usr/local/bin/before-notebook.d/

# TODO: Remove PGSQL and daemon log files, and other unneeded files
RUN --mount=from=qe_conda_env,source=${QE_DIR},target=${QE_DIR} \
    bash /usr/local/bin/before-notebook.d/20_start-postgresql.sh && \
    bash /usr/local/bin/before-notebook.d/40_prepare-aiida.sh && \
    bash /usr/local/bin/before-notebook.d/42_setup-hq-computer.sh && \
    python -m aiidalab_qe install-qe --computer ${COMPUTER_LABEL} && \
    python -m aiidalab_qe install-pseudos --source ${PSEUDO_FOLDER} && \
    # setup code: pythonjob, bader, wannier90 code
    verdi code create core.code.installed --label python --computer=${COMPUTER_LABEL} --default-calc-job-plugin pythonjob.pythonjob --filepath-executable=/opt/conda/bin/python -n && \
    verdi code create core.code.installed --label bader --computer=${COMPUTER_LABEL} --default-calc-job-plugin bader.bader --filepath-executable=${QE_DIR}/bin/bader -n && \
    verdi code create core.code.installed --label wannier90 --computer=${COMPUTER_LABEL} --default-calc-job-plugin wannier90.wannier90 --filepath-executable=${QE_DIR}/bin/wannier90.x -n && \
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
    cd /home/${NB_USER} && tar -cf /opt/conda/home.tar \
      --exclude .npm/_logs --exclude .aiida/access/default --exclude .cache --exclude .npm/ --exclude work --exclude .conda .

###############################################################################
# 7) Prepare a clean copy of the QeApp repo
#    - This must be a separate stage because we run `git clean` after COPY,
#      and we don't want to have junk in intermediate layers.
###############################################################################
FROM base AS qeapp_src

WORKDIR ${QE_APP_FOLDER}
COPY --chown=${NB_UID}:${NB_GID} . .
# Remove all untracked files and directories.
RUN git clean -dffx

###############################################################################
# 8) Final stage
#    - Installs python dependencies again, copies QE env, compiles wannier90,
#      and sets up the final environment.
###############################################################################
FROM base
ARG QE_DIR
ARG QE_APP_SRC
ARG TARGETARCH

USER root
# Install deps for wannier90
RUN apt-get -q update && \
    apt-get -q install -y --no-install-recommends liblapack3 openmpi-bin && \
    apt-get -q clean && \
    rm -rf /var/lib/apt/lists/*

USER ${NB_USER}
WORKDIR /tmp

# Install common dependencies such as pymatgen and pandas via conda,
# to avoid building them when installing via pip
# TODO: Remove this once it is part of the full-stack image.
RUN mamba install aiida-core.atomic_tools -y && \
    mamba clean --all -f -y

RUN python -m pip install --no-user --no-cache-dir ${AIIDA_HQ_PKG}

# copy hq binary
COPY --from=home_build /opt/conda/hq /usr/local/bin/
# Copy the QE conda environment
COPY --from=qe_conda_env ${QE_DIR} ${QE_DIR}
# This is for backwards compatibility
COPY --from=qe_conda_env ${QE_DIR}/bin/wannier90.x /opt/conda/bin/wannier90.x

USER root

# We exclude 42_setup-hq-computer.sh file because the computer is already setup, thus it is not needed in the final image.
COPY ./before-notebook.d/00_untar-home.sh ./before-notebook.d/43_start-hq.sh /usr/local/bin/before-notebook.d/

# Remove the content of /home/<user>, but keep the folder itself
RUN find /home/${NB_USER}/ -mindepth 1 -delete

ENV HOME_TAR="/opt/home.tar"
COPY --from=home_build /opt/conda/home.tar "$HOME_TAR"

COPY --from=qeapp_src ${QE_APP_FOLDER} ${QE_APP_FOLDER}

USER ${NB_USER}
WORKDIR "/home/${NB_USER}"
