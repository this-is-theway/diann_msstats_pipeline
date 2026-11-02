# Minimal ThermoRawFileParser image for the RAW_TO_MZML process.
#
# Build:  docker build -t yourregistry/thermorawfileparser:1.4.5 .

FROM mambaorg/micromamba:1.5.10

USER root

# Pin explicitly — check https://anaconda.org/bioconda/thermorawfileparser
# for the current version list before bumping this.
ARG TRFP_VERSION=1.4.5

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        thermorawfileparser=${TRFP_VERSION} \
    && micromamba clean --all --yes


ENV PATH=/opt/conda/bin:$PATH
ENTRYPOINT []