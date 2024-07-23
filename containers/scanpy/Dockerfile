# syntax=docker/dockerfile:1

ARG PYTHON_VERSION=3.11.9-slim-bookworm
FROM python:${PYTHON_VERSION} AS base

# Prevents Python from writing pyc files.
ENV PYTHONDONTWRITEBYTECODE=1

# Keeps Python from buffering stdout and stderr to avoid situations where
# the application crashes without emitting any logs due to buffering.
ENV PYTHONUNBUFFERED=1

LABEL python_version=${PYTHON_VERSION}
LABEL authors="Michael Weinberger <michael.weinberger11@gmail.com>"

# install packages
RUN apt-get update && \
	apt-get install -y procps && \
    apt-get clean -y && \
    rm -rf /var/lib/apt/lists/*

# Declare python packages to be installed
ENV PYTHON_PACKAGES="\
    anndata==0.9.1 \
    doubletdetection==4.2 \
    matplotlib==3.7.1 \
    numpy==1.24.3 \
    pandas==2.1.4 \
    scanpy==1.9.3 \
    scipy==1.10.1 \
    seaborn==0.12.2 \
"

# Download dependencies as a separate step to take advantage of Docker's caching.
# Leverage a bind mount to requirements.txt to avoid having to copy them into
# into this layer.
RUN python -m pip install --no-cache-dir ${PYTHON_PACKAGES}

# Done due to an error introduced by new version of scipy
#RUN python -m pip install --no-cache-dir --force-reinstall --no-deps "scipy<1.14.0"

# Create a non-privileged user that the app will run under.
# See https://docs.docker.com/go/dockerfile-user-best-practices/
ARG UID=10001
RUN adduser \
    --disabled-password \
    --gecos "" \
    --home "/nonexistent" \
    --shell "/sbin/nologin" \
    --no-create-home \
    --uid "${UID}" \
    appuser

# Switch to the non-privileged user to run the application.
USER appuser

# by default /bin/bash is executed
#ENTRYPOINT ["/bin/bash"]