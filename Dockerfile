FROM nvidia/cuda:12.2.0-devel-ubuntu22.04

SHELL ["/bin/bash", "-c"]

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# 1) System packages
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
    libpq-dev \
    cmake \
    nano \
    postgresql-client \
    && rm -rf /var/lib/apt/lists/*

# 2) Add CRAN repo for R 4.3.x
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | apt-key add - \
 && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/"

RUN apt-get update && apt-get install -y --no-install-recommends r-base \
    && rm -rf /var/lib/apt/lists/*

# 3) Install Miniconda
ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
 && /bin/bash /tmp/miniconda.sh -b -p $CONDA_DIR \
 && rm /tmp/miniconda.sh
ENV PATH=$CONDA_DIR/bin:$PATH

# 4) Create conda environment from environment.yml
COPY environment.yml /tmp/environment.yml
RUN conda update -n base conda -y \
 && conda env create -f /tmp/environment.yml -n myenv \
 && conda clean -afy

# Make sure conda env is on PATH
ENV PATH="/opt/conda/envs/myenv/bin:$PATH"

# 5) Install R "renv" + copy your entire code into /app
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org')" \
 && R -e "remotes::install_github('rstudio/renv')"

WORKDIR /app
COPY . /app

# 6) renv::restore() in /app
RUN R -e "setwd('/app'); renv::restore(confirm = FALSE, library = '/app/renv/library')"

# 7) Optional: Build GROMACS
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

# 8) Expose port if needed for Streamlit
EXPOSE 8501

CMD ["/bin/bash"]
