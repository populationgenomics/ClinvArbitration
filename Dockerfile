# Any analysis-runner driver image must at least include git.
ARG PY_VER=${PY_VER:-3.10}

FROM python:${PY_VER}-slim-bullseye AS basic

ENV PYTHONDONTWRITEBYTECODE=1

# ClinvArbitration pipeline version.
ENV VERSION=2.2.6

RUN apt update && apt install --no-install-recommends -y \
        apt-transport-https \
        bzip2 \
        ca-certificates \
        git \
        gnupg \
        libbz2-1.0 \
        libcurl4 \
        liblzma5 \
        openjdk-17-jdk-headless \
        wget \
        zip && \
    pip install --no-cache-dir -U pip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

FROM basic AS bcftools_compiler

RUN apt-get update && apt-get install --no-install-recommends -y \
        gcc \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libssl-dev \
        make \
        zlib1g-dev && \
    wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 && \
    tar -xf bcftools-1.22.tar.bz2 && \
    cd bcftools-1.22 && \
    ./configure --enable-libcurl --enable-s3 --enable-gcs && \
    make && \
    strip bcftools plugins/*.so && \
    make DESTDIR=/bcftools_install install

FROM basic AS base_bcftools

COPY --from=bcftools_compiler /bcftools_install/usr/local/bin/* /usr/local/bin/
COPY --from=bcftools_compiler /bcftools_install/usr/local/libexec/bcftools/* /usr/local/libexec/bcftools/

FROM base_bcftools AS now_build_clinvarbitration

# now do some fun stuff, installing ClinvArbitration
WORKDIR /clinvarbitration

COPY src src/
COPY LICENSE pyproject.toml README.md ./

# pip install but don't retain the cache files
RUN pip install --no-cache-dir .

COPY nextflow nextflow/
