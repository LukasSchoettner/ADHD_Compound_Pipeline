FROM nvidia/cuda:12.6.3-devel-ubuntu24.04

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
    gcc-13 \
    g++-13 \
    ca-certificates \
    gfortran \
    libblas-dev \
    liblapack-dev \
    libssl-dev \
    libxml2-dev \
    libpq-dev \
    libtool \
    autoconf \
    automake \
    cmake \
    zlib1g-dev \
    nano \
    postgresql-client \
    && rm -rf /var/lib/apt/lists/*

# 2) Install Miniconda
ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
 && /bin/bash /tmp/miniconda.sh -b -p $CONDA_DIR \
 && rm /tmp/miniconda.sh
ENV PATH=$CONDA_DIR/bin:$PATH

# 3) Create conda environment from environment.yml
COPY environment.yml /tmp/environment.yml
RUN conda update -n base conda -y \
 && conda env create -f /tmp/environment.yml -n myenv \
 && conda clean -afy

# Make sure Conda environment is available in the shell
#ENV PATH="/opt/conda/envs/myenv/bin:$CONDA_DIR/bin:$PATH"
RUN echo "source /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc \
 && echo "conda activate myenv" >> ~/.bashrc

# Make sure conda env is on PATH
#ENV PATH="/opt/conda/envs/myenv/bin:$PATH"
#ENV LD_LIBRARY_PATH="/opt/conda/envs/myenv/lib:$LD_LIBRARY_PATH"

ENV PATH="/opt/conda/envs/myenv/bin:$CONDA_DIR/bin:$PATH"
ENV LD_LIBRARY_PATH="/opt/conda/envs/myenv/lib:$LD_LIBRARY_PATH"
ENV PKG_CONFIG_PATH="/opt/conda/envs/myenv/lib/pkgconfig:$PKG_CONFIG_PATH"

# 4) Install R and necessary packages via Conda
#RUN conda install -n myenv -c conda-forge \
#    r-base \
#    r-curl \
#    r-renv \
#    && conda clean -afy

# 5) Install R "renv" + copy your entire code into /app
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org')" \
 && R -e "remotes::install_github('rstudio/renv')"

# 6) Copy your entire code into /app
WORKDIR /app
COPY . /app

# 7) renv::restore() in /app
RUN R -e "renv::restore(confirm = FALSE, library = '/app/renv/library')"
# Add before the GROMACS build step
ENV CC=/usr/bin/gcc-13
ENV CXX=/usr/bin/g++-13

# 8) Optional: Build GROMACS
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

# 9) Expose port if needed for Streamlit
EXPOSE 8501

CMD ["bash", "-l"]
