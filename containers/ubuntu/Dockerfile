# syntax=docker/dockerfile:1

ARG UBUNTU_VER=22.04
FROM ubuntu:${UBUNTU_VER}

LABEL authors="Michael Weinberger <michael.weinberger11@gmail.com>"

# install packages
RUN apt-get update && \
	apt-get install -y build-essential wget gzip rsync procps && \
    apt-get clean -y && \
    rm -rf /var/lib/apt/lists/*

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