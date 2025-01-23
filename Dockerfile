# -----------------------------------------------------------------------------
# 1) Use an NVIDIA CUDA base image (GPU example) or a Miniconda base (CPU only)
# -----------------------------------------------------------------------------
FROM nvidia/cuda:12.2.0-devel-ubuntu22.04

SHELL ["/bin/bash", "-c"]

# -----------------------------------------------------------------------------
# 2) Set environment variables to avoid tzdata interactive prompt
# -----------------------------------------------------------------------------
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# -----------------------------------------------------------------------------
# 3) Prepare system to install latest R from CRAN
#    (so we get R >= 4.3 instead of Ubuntu's older R 4.1.2)
# -----------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
    gnupg2 \
    wget \
    git \
    build-essential \
    ca-certificates \
    gfortran \
    libblas-dev \
    liblapack-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Add CRAN repository for Ubuntu jammy (22.04)
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | apt-key add - \
 && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/"

RUN apt-get update && apt-get install -y --no-install-recommends r-base \
    && rm -rf /var/lib/apt/lists/*

# -----------------------------------------------------------------------------
# 4) Install Miniconda
# -----------------------------------------------------------------------------
ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
 && /bin/bash /tmp/miniconda.sh -b -p $CONDA_DIR \
 && rm /tmp/miniconda.sh
ENV PATH=$CONDA_DIR/bin:$PATH

# -----------------------------------------------------------------------------
# 5) Copy and create your conda environment from environment.yml
# -----------------------------------------------------------------------------
COPY environment.yml /tmp/environment.yml

RUN conda update -n base conda -y \
 && conda env create -f /tmp/environment.yml -n myenv \
 && conda clean -afy

# Activate 'myenv' by default
RUN echo "conda activate myenv" >> ~/.bashrc

# -----------------------------------------------------------------------------
# 6) R + renv setup
# -----------------------------------------------------------------------------
COPY renv.lock /tmp/renv.lock

# Install "remotes" & "renv" into R 4.3.x
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org')" \
 && R -e "remotes::install_github('rstudio/renv')"

# Copy your entire repo (including renv/ if present) into /workspace
WORKDIR /workspace
COPY . /workspace

# Run renv::restore -- now we have R 4.3.x, so Bioc 3.18 should work
RUN R -e "setwd('/workspace'); renv::restore(confirm = FALSE)"

# -----------------------------------------------------------------------------
# 7) Build GROMACS with CUDA + RDTSCP
# -----------------------------------------------------------------------------
ARG GROMACS_VERSION="v2024.4"
RUN git clone --branch ${GROMACS_VERSION} https://gitlab.com/gromacs/gromacs.git /tmp/gromacs \
 && mkdir -p /tmp/gromacs/build \
 && cd /tmp/gromacs/build \
 && cmake .. \
      -DGMX_CUDA=ON \
      -DGMX_USE_RDTSCP=ON \
      -DGMX_BUILD_OWN_FFTW=ON \
      -DGMX_MPI=OFF \
      -DGMX_OPENMP=ON \
      -DCMAKE_INSTALL_PREFIX=/opt/gromacs \
      -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda \
 && make -j 8 \
 && make install \
 && rm -rf /tmp/gromacs

ENV PATH="/opt/gromacs/bin:$PATH"

# -----------------------------------------------------------------------------
# 8) Expose Streamlit or Shiny ports, if needed
# -----------------------------------------------------------------------------
EXPOSE 8501

# -----------------------------------------------------------------------------
# 9) Default command
# -----------------------------------------------------------------------------
CMD ["/bin/bash"]