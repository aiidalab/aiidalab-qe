# ##
# # Unified multi-arch Dockerfile for both AMD64 and ARM64 builds.
# # Conditionally installs 'bader' in conda (on x86_64) or compiles from source (on arm64),
# # and downloads the correct Hyperqueue (hq) binary for each architecture.
# ##

# ###############################################################################
# # 1) Global ARGs
# ###############################################################################
ARG FULL_STACK_VER=2025.1026
# ARG UV_VER=0.4.7
# ARG QE_VER=7.4
# ARG QE_DIR=/opt/conda/envs/quantum-espresso-${QE_VER}
# ARG HQ_VER=0.19.0

# ARG UV_CACHE_DIR=/tmp/uv_cache
# ARG QE_APP_SRC=/tmp/quantum-espresso
# ARG COMPUTER_LABEL="localhost"

# #
# # We'll define the possible HQ download URLs (for x86_64 and ARM64).
# #
# ARG HQ_URL_AMD64="https://github.com/It4innovations/hyperqueue/releases/download/v${HQ_VER}/hq-v${HQ_VER}-linux-x64.tar.gz"
# ARG HQ_URL_ARM64="https://github.com/It4innovations/hyperqueue/releases/download/v${HQ_VER}/hq-v${HQ_VER}-linux-arm64-linux.tar.gz"

# ###############################################################################
# # 2) uv stage (unchanged)
# ###############################################################################
# FROM ghcr.io/astral-sh/uv:${UV_VER} AS uv

# ###############################################################################
# # 3) qe_conda_env stage
# #    - Creates the quantum-espresso conda environment in /opt/conda/envs/quantum-espresso-<QE_VER>
# #    - On x86_64, we install the 'bader' package from conda. On arm64, skip 'bader' because
# #      it’s unavailable and we'll build from source later.
# ###############################################################################
# FROM ghcr.io/aiidalab/full-stack:${FULL_STACK_VER} AS qe_conda_env
# ARG QE_VER
# ARG QE_DIR
# # Docker sets TARGETARCH automatically (e.g. "amd64" or "arm64")
# ARG TARGETARCH

# USER ${NB_USER}
# RUN set -ex; \
#     if [ "${TARGETARCH}" = "amd64" ]; then \
#       echo "Installing QE plus Bader for x86_64..."; \
#       mamba create -p ${QE_DIR} --yes qe=${QE_VER} bader; \
#     else \
#       echo "Installing QE (without bader) for ARM64..."; \
#       mamba create -p ${QE_DIR} --yes qe=${QE_VER}; \
#     fi && \
#     mamba clean --all -f -y

# ###############################################################################
# # 4) build_deps stage
# #    - Installs Python dependencies using uv for caching
# ###############################################################################
# FROM ghcr.io/aiidalab/full-stack:${FULL_STACK_VER} AS build_deps
# ARG QE_DIR
# ARG UV_CACHE_DIR
# ARG QE_APP_SRC

# WORKDIR ${QE_APP_SRC}

# COPY --chown=${NB_UID}:${NB_GID} src/ ${QE_APP_SRC}/src
# COPY --chown=${NB_UID}:${NB_GID} setup.cfg pyproject.toml LICENSE README.md ${QE_APP_SRC}

# ENV UV_CONSTRAINT=${PIP_CONSTRAINT}
# RUN --mount=from=uv,source=/uv,target=/bin/uv \
#     uv pip install --strict --system --cache-dir=${UV_CACHE_DIR} .

# ###############################################################################
# # 5) home_build stage
# #    - Prepares AiiDA profile, sets up hyperqueue, installs QE codes/pseudos,
# #      and archives the home folder (home.tar).
# ###############################################################################
# FROM build_deps AS home_build
# ARG UV_CACHE_DIR
# ARG QE_DIR
# ARG HQ_VER
# ARG COMPUTER_LABEL
# # We'll use these to pick the correct HQ binary
# ARG TARGETARCH
# ARG HQ_URL_AMD64
# ARG HQ_URL_ARM64

# #
# # Download and unpack the correct hq binary for the architecture:
# #
# RUN set -ex; \
#     if [ "${TARGETARCH}" = "arm64" ]; then \
#       echo "Downloading hyperqueue for ARM64..."; \
#       wget -c -O hq.tar.gz "${HQ_URL_ARM64}"; \
#     else \
#       echo "Downloading hyperqueue for x86_64..."; \
#       wget -c -O hq.tar.gz "${HQ_URL_AMD64}"; \
#     fi && \
#     tar xf hq.tar.gz -C /opt/conda/

# ENV PSEUDO_FOLDER=/tmp/pseudo

# RUN mkdir -p ${PSEUDO_FOLDER} && \
#     python -m aiidalab_qe download-pseudos --dest ${PSEUDO_FOLDER}

# ENV UV_CONSTRAINT=${PIP_CONSTRAINT}
# # Install the aiida-hyperqueue
# # XXX: fix me after release aiida-hyperqueue
# RUN --mount=from=uv,source=/uv,target=/bin/uv \
#     --mount=from=build_deps,source=${UV_CACHE_DIR},target=${UV_CACHE_DIR},rw \
#     uv pip install --system --strict --cache-dir=${UV_CACHE_DIR} \
#     "aiida-hyperqueue@git+https://github.com/aiidateam/aiida-hyperqueue"

# COPY ./before-notebook.d/* /usr/local/bin/before-notebook.d/

# # TODO: Remove PGSQL and daemon log files, and other unneeded files
# RUN --mount=from=qe_conda_env,source=${QE_DIR},target=${QE_DIR} \
#     bash /usr/local/bin/before-notebook.d/20_start-postgresql.sh && \
#     bash /usr/local/bin/before-notebook.d/40_prepare-aiida.sh && \
#     bash /usr/local/bin/before-notebook.d/42_setup-hq-computer.sh && \
#     python -m aiidalab_qe install-qe --computer ${COMPUTER_LABEL} && \
#     python -m aiidalab_qe install-pseudos --source ${PSEUDO_FOLDER} && \
#     # steup code: pythonjob, bader, wannier90 code
#     verdi code create core.code.installed --label python --computer=localhost --default-calc-job-plugin pythonjob.pythonjob --filepath-executable=/opt/conda/bin/python -n && \
#     verdi code create core.code.installed --label bader --computer=localhost --default-calc-job-plugin bader.bader --filepath-executable=${QE_DIR}/bin/bader -n && \
#     verdi code create core.code.installed --label wannier90 --computer=localhost --default-calc-job-plugin wannier90.wannier90 --filepath-executable=/opt/conda/bin/wannier90.x -n && \
#     # Additional plugins
#     pip uninstall -y phonopy && \
#     pip install aiida-bader \
#     git+https://github.com/mikibonacci/aiidalab-qe-vibroscopy@v1.2.0 \
#     git+https://github.com/mikibonacci/aiidalab-qe-muon@v1.0.0 && \
#     # run post_install for plugin
#     python -m aiida_bader post-install && \
#     python -m aiidalab_qe_vibroscopy setup-phonopy && \
#     python -m aiidalab_qe_muon setup-python3 && \
#     # wannier90 plugin need SSSP 1.1
#     aiida-pseudo install sssp -v 1.1 -x PBE && \
#     aiida-pseudo install sssp -v 1.1 -x PBEsol && \
#     verdi daemon stop && \
#     mamba run -n aiida-core-services pg_ctl stop && \
#     touch /home/${NB_USER}/.FLAG_HOME_INITIALIZED && \
#     # NOTE: The work folder is empty but if included clashes with the work folder in a Renku
#     # session whose permissions cannot be changed.
#     # For the same permisssion reason, the .conda folder clashes with aiidalab-launch.
#     # It is usually safe (and preferable) to let .conda be recreated on the fly each time,
#     # because .conda typically just holds local environment information, caches, or references
#     # to available environments.
#     cd /home/${NB_USER} && tar -cf /opt/conda/home.tar --exclude work --exclude .conda .

###############################################################################
# 6) Final stage
#    - Installs python dependencies again, copies QE env, compiles wannier90,
#      (conditionally) compiles bader on ARM64, and sets up the final environment.
###############################################################################
FROM ghcr.io/aiidalab/full-stack:${FULL_STACK_VER}
ARG QE_DIR
ARG QE_APP_SRC
ARG UV_CACHE_DIR
ARG COMPUTER_LABEL
ARG TARGETARCH

# USER ${NB_USER}
# WORKDIR /tmp
# # Install python dependencies
# # Use uv cache from the previous build step
# # # Install the aiida-hyperqueue
# # # XXX: fix me after release aiida-hyperqueue
# ENV UV_CONSTRAINT=${PIP_CONSTRAINT}
# RUN --mount=from=uv,source=/uv,target=/bin/uv \
#     --mount=from=build_deps,source=${UV_CACHE_DIR},target=${UV_CACHE_DIR},rw \
#     --mount=from=build_deps,source=${QE_APP_SRC},target=${QE_APP_SRC},rw \
#     uv pip install --strict --system --compile-bytecode --cache-dir=${UV_CACHE_DIR} ${QE_APP_SRC} "aiida-hyperqueue@git+https://github.com/aiidateam/aiida-hyperqueue"
# # Install plugins in the final image
# RUN pip install aiida-bader \
#     git+https://github.com/mikibonacci/aiidalab-qe-vibroscopy@v1.2.0 \
#     git+https://github.com/mikibonacci/aiidalab-qe-muon@v1.0.0

# # copy hq binary
# COPY --from=home_build /opt/conda/hq /usr/local/bin/
# # Copy the QE conda environment
# COPY --from=qe_conda_env ${QE_DIR} ${QE_DIR}

# USER root
# RUN apt-get update && \
#     apt-get install -y --no-install-recommends \
#       gfortran libblas-dev liblapack-dev git openmpi-bin

# # Build wannier90 for all arches, and build bader from source ONLY on arm64
# RUN set -ex; \
#     git clone --depth=1 https://github.com/wannier-developers/wannier90.git /tmp/wannier90; \
#     cd /tmp/wannier90; \
#     cp config/make.inc.gfort make.inc; \
#     make wannier; \
#     cp wannier90.x /opt/conda/bin/wannier90.x; \
#     \
#     if [ "${TARGETARCH}" = "arm64" ]; then \
#       echo "Building bader from source for ARM64..."; \
#       git clone --depth=1 https://gitlab.com/jameskermode/bader.git /tmp/bader; \
#       cd /tmp/bader; \
#       cp makefile.osx_gfortran Makefile; \
#       make; \
#       cp bader ${QE_DIR}/bin/bader; \
#     else \
#       echo "Skipping Bader build on AMD64 (installed via conda)."; \
#     fi; \
#     \
#     apt-get remove --purge -y gfortran libblas-dev liblapack-dev && \
#     apt-get install -y --no-install-recommends liblapack3 && \
#     apt-get autoremove -y && \
#     apt-get clean && \
#     rm -rf /var/lib/apt/lists/* /tmp/wannier90 /tmp/bader

# # We exclude 42_setup-hq-computer.sh file because the computer is already steup, thus it is not needed in the final image.
# COPY ./before-notebook.d/00_untar-home.sh ./before-notebook.d/43_start-hq.sh /usr/local/bin/before-notebook.d/
# ENV COMPUTER_LABEL=$COMPUTER_LABEL

# # Remove the content of /home/<user>, but keep the folder itself
# RUN find /home/${NB_USER}/ -mindepth 1 -delete

# ENV QE_APP_FOLDER=/opt/conda/quantum-espresso
# COPY --chown=${NB_UID}:${NB_GID} . ${QE_APP_FOLDER}
# # Remove all untracked files and directories.
# RUN git clean -dffx || true

# ENV HOME_TAR="/opt/home.tar"
# COPY --from=home_build /opt/conda/home.tar "$HOME_TAR"

USER ${NB_USER}
WORKDIR "/home/${NB_USER}"
