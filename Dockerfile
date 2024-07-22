# syntax=docker/dockerfile:1
FROM ghcr.io/astral-sh/uv:0.2.18 AS uv

FROM ghcr.io/aiidalab/full-stack:2024.1021 AS home_build
USER ${NB_USER}

ENV QE_APP_FOLDER=/opt/conda/quantum-espresso
WORKDIR ${QE_APP_FOLDER}

COPY --chown=${NB_UID}:${NB_GID} src/ ${QE_APP_FOLDER}/src
COPY --chown=${NB_UID}:${NB_GID} setup.cfg pyproject.toml *yaml README.md ${QE_APP_FOLDER}

ENV UV_CONSTRAINT=${PIP_CONSTRAINT}
RUN --mount=from=uv,source=/uv,target=/bin/uv \
    uv pip install --system --no-cache . && \
    rm -rf build/ src/aiidalab_qe.egg-info/

# 4. Prepare AiiDA profile and localhost computer
# 5. Install the QE pseudopotentials and codes
# TODO: Remove PGSQL and daemon log files, and other unneeded files
RUN bash /usr/local/bin/before-notebook.d/20_start-postgresql.sh && \
    bash /usr/local/bin/before-notebook.d/40_prepare-aiida.sh && \
    python -m aiidalab_qe install-qe && \
    python -m aiidalab_qe install-pseudos && \
    verdi daemon stop && \
    mamba run -n aiida-core-services pg_ctl stop && \
    cd /home/${NB_USER} && tar -cf /opt/conda/home.tar .

FROM ghcr.io/aiidalab/full-stack:2024.1021
COPY --from=home_build /opt/conda/home.tar /opt/conda

USER ${NB_USER}

ENV QE_VERSION="7.2"
# 1. Install Quantum Espresso into a conda environment
# (we do this here for better caching during Docker build)
RUN mamba create -p /opt/conda/envs/quantum-espresso --yes \
    qe=${QE_VERSION} && \
    mamba clean --all -f -y

ENV QE_APP_FOLDER=/opt/conda/quantum-espresso
WORKDIR "${QE_APP_FOLDER}"

COPY --chown=${NB_UID}:${NB_GID} . ${QE_APP_FOLDER}
# Remove all untracked files and directories.
RUN git clean -dffx || true

# 3.Install python dependencies
# Use uv instead of pip to speed up installation, per docs:
# https://github.com/astral-sh/uv/blob/main/docs/guides/docker.md#using-uv-temporarily
# Use the same constraint file as pip
ENV UV_CONSTRAINT=${PIP_CONSTRAINT}
RUN --mount=from=uv,source=/uv,target=/bin/uv \
    uv pip install --system --no-cache . && \
    rm -rf build/ src/aiidalab_qe.egg-info/

USER root
COPY ./before-notebook.d/* /usr/local/bin/before-notebook.d/
RUN fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

WORKDIR "/home/${NB_USER}"
USER ${NB_USER}
