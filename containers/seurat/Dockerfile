# syntax=docker/dockerfile:1

ARG SEURAT_VERSION=4.3.0
FROM satijalab/seurat:${SEURAT_VERSION}

LABEL authors="Michael Weinberger <michael.weinberger11@gmail.com>"

# install packages
RUN apt-get update && \
	apt-get install -y libcurl4-openssl-dev libssl-dev procps && \
    apt-get clean -y && \
    rm -rf /var/lib/apt/lists/*

# Install R dependencies
RUN install2.r --error \ 
    dplyr \
    ggplot2 \
    patchwork \
    pheatmap \
    remotes \
    viridis \
    writexl \
&& R -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')" \
&& R -e "remotes::install_github('immunogenomics/harmony')"

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

# Switch to the non-privileged user
USER appuser

# by default /bin/bash is executed
#ENTRYPOINT ["/bin/bash"]