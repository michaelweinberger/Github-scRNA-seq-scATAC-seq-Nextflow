# syntax=docker/dockerfile:1

ARG UBUNTU_VER=22.04
FROM ubuntu:${UBUNTU_VER}

LABEL authors="Michael Weinberger <michael.weinberger11@gmail.com>"

WORKDIR /

# install packages
RUN apt-get update && \
	apt-get install -y wget vim bzip2 gcc procps && \
    apt-get clean -y && \
    rm -rf /var/lib/apt/lists/*

# install miniconda + samtools
ARG MINICONDA_VER=py311_24.5.0-0
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA_VER}-Linux-x86_64.sh -O ~/miniconda.sh && \
	bash ~/miniconda.sh -b -p /opt/conda && \
	rm ~/miniconda.sh

ENV PATH=/opt/conda/bin:$PATH

RUN conda update conda -y && \
    conda install libgcc-ng --channel conda-forge && \
    conda install samtools=1.17 --channel bioconda

# Prevents Python from writing pyc files.
ENV PYTHONDONTWRITEBYTECODE=1

# Keeps Python from buffering stdout and stderr to avoid situations where
# the application crashes without emitting any logs due to buffering.
ENV PYTHONUNBUFFERED=1

# Install further python packages
ENV PYTHON_PACKAGES="\
    click==8.1.3 \
    cython==0.29.34 \
    h5py==3.8.0 \
    matplotlib==3.7.1 \
    numba==0.58.1 \
    scikit-learn==1.4.2 \
    scipy==1.10.1 \
"
RUN python -m pip install --no-cache-dir ${PYTHON_PACKAGES} && \
    python -m pip install --no-cache-dir velocyto==0.17.17

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