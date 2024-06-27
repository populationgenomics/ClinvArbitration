FROM python:3.10-bullseye

# take as a command line argument, or
ARG RELEASE=${RELEASE:-1.2.0}

RUN apt update && apt install -y \
        apt-transport-https \
        bzip2 \
        ca-certificates \
        git \
        gnupg \
        openjdk-11-jdk-headless \
        wget \
        zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

RUN pip install --no-cache-dir git+https://github.com/populationgenomics/ClinvArbitration.git@${RELEASE}
