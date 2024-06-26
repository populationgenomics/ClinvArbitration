FROM python:3.10-bullseye

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

COPY requirements*.txt .

RUN pip install -r requirements.txt
COPY README.md .
COPY setup.py .
COPY clinvarbitration clinvarbitration/

RUN pip install --no-cache-dir .
