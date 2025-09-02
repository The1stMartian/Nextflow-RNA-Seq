FROM mambaorg/micromamba:1.5.8
ARG MAMBA_DOCKERFILE_ACTIVATE=1
SHELL ["/bin/bash", "-lc"]

ENV DEBIAN_FRONTEND=noninteractive LC_ALL=C.UTF-8 LANG=C.UTF-8 TZ=UTC

USER root

# System deps, Java for Nextflow, and bio tools from apt
RUN apt-get update && apt-get install -y --no-install-recommends \
      locales bash curl wget git unzip bzip2 make g++ \
      libxml2-dev libssl-dev libcurl4-openssl-dev \
      libfontconfig1-dev libfreetype6-dev libpng-dev libtiff-dev libjpeg-dev \
      libharfbuzz-dev libfribidi-dev libglpk-dev libgit2-dev \
      openjdk-17-jre-headless \
      # bio tools:
      samtools subread rna-star \
      # python (kept for scripts, but no pip build needed for umi-tools now)
      python3 python3-pip python3-setuptools python3-wheel \
  && echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
  && locale-gen \
  && rm -rf /var/lib/apt/lists/*


# Install minimal system packages
RUN apt-get update && apt-get install -y --no-install-recommends \
      curl ca-certificates unzip procps pigz && \
    rm -rf /var/lib/apt/lists/*
	
# R and packages (CRAN + Bioconductor)
RUN apt-get update && apt-get install -y --no-install-recommends \
        software-properties-common \
        dirmngr \
        gnupg \
        build-essential \
        libcurl4-gnutls-dev \
        libxml2-dev \
        libssl-dev \
        libssh2-1-dev \
        pandoc 
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
        add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
        apt-get update && \
        apt-get install -y --no-install-recommends \
        r-base \
        r-base-dev \
        r-recommended
RUN R -e "install.packages(c('readr','dplyr','tibble','ggplot2','stringr','tidyverse','ggrepel','scales','BiocManager','purrr'), repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install('DESeq2', ask=FALSE, update=TRUE)"

# Configure micromamba channels
RUN micromamba config set channel_priority strict && \
    micromamba config append channels conda-forge && \
    micromamba config append channels bioconda && \
    micromamba config append channels defaults

# Get UMI-tools
RUN micromamba install umi_tools

# Install AWS CLI v2 (outside conda)
RUN curl -fsSL "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "/tmp/awscliv2.zip" && \
    unzip -q /tmp/awscliv2.zip -d /tmp && \
    /tmp/aws/install && rm -rf /tmp/aws /tmp/awscliv2.zip

# Add conda env to PATH for Nextflow processes
ENV PATH="/opt/conda/envs/bio/bin:/usr/local/bin:${PATH}"

# Symlink Java for tools that expect it in /usr/bin
RUN ln -sf /opt/conda/envs/bio/bin/java /usr/bin/java || true

# Create non-root user for Nextflow
ARG USERNAME=nfuser UID=1000
RUN useradd -m -u ${UID} -s /bin/bash ${USERNAME} && \
    mkdir -p /workspace && chown ${USERNAME}:${USERNAME} /workspace
USER ${USERNAME}
WORKDIR /workspace

# No ENTRYPOINT needed; Nextflow will run commands directly