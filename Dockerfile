# syntax=docker/dockerfile:1
ARG FULL_STACK_VER=2025.1026
ARG UV_VER=0.4.7
ARG QE_VER=7.4
ARG QE_DIR=/opt/conda/envs/quantum-espresso-${QE_VER}
ARG HQ_VER=0.19.0

FROM ghcr.io/aiidalab/full-stack:${FULL_STACK_VER}
ARG QE_DIR
ARG QE_APP_SRC
ARG UV_CACHE_DIR
ARG COMPUTER_LABEL
USER ${NB_USER}

WORKDIR "/home/${NB_USER}"
