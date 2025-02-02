FROM nvidia/cuda:12.3.1-devel-ubuntu22.04

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
    gcc-12 \
    g++-12 \
    ca-certificates \
    gfortran \
    libblas-dev \
    liblapack-dev \
    libssl-dev \
    libxml2-dev \
    libpq-dev \
    libtool \
    libgl1 \
    autoconf \
    automake \
    cmake \
    zlib1g-dev \
    nano \
    postgresql-client \
    swig \
    make \
    libeigen3-dev \
    libopenbabel-dev \
    libncurses6 \
    libtinfo6 \
    && rm -rf /var/lib/apt/lists/*

# Set GCC 12 and G++ 12 as default
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-12 100 \
    && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-12 100 \
    && update-alternatives --set gcc /usr/bin/gcc-12 \
    && update-alternatives --set g++ /usr/bin/g++-12

# Create /lib64 symlink as a workaround for linker errors
RUN ln -s /lib/x86_64-linux-gnu /lib64

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

# Ensure Conda environment is activated in the shell
RUN echo "source /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc \
 && echo "conda activate myenv" >> ~/.bashrc

ENV PATH="/opt/conda/envs/myenv/bin:$PATH"
# Removed LD_LIBRARY_PATH to prevent library conflicts
# ENV PKG_CONFIG_PATH="/opt/conda/envs/myenv/lib/pkgconfig:$PKG_CONFIG_PATH"

# 4) Verify Open Babel installation
RUN obabel --version || true

# 5) Copy project code into /app
WORKDIR /app
COPY . /app

# 6) Restore R environment (if applicable)
RUN R -e "renv::restore(confirm = FALSE, library = '/app/renv/library')"

# 7) Add GROMACS (Optional)
ARG GROMACS_VERSION="v2024.4"
RUN git clone --branch ${GROMACS_VERSION} https://gitlab.com/gromacs/gromacs.git /tmp/gromacs \
    && mkdir -p /tmp/gromacs/build \
    && cd /tmp/gromacs/build \
    && cmake .. \
        -DGMX_GPU=CUDA \
        -DGMX_USE_RDTSCP=ON \
        -DFORCE_EXTERNAL_FFTW=ON \
        -DGMX_MPI=OFF \
        -DGMX_OPENMP=ON \
        -DCMAKE_INSTALL_PREFIX=/opt/gromacs \
        -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda \
        -DCMAKE_PREFIX_PATH=/opt/conda/envs/myenv \
        -DCMAKE_C_COMPILER=/usr/bin/gcc-12 \
        -DCMAKE_CXX_COMPILER=/usr/bin/g++-12 \
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
        -DCMAKE_CUDA_HOST_COMPILER=/usr/bin/gcc-12 \
    && make -j$(nproc) VERBOSE=1 \
    && make install \
    && rm -rf /tmp/gromacs

ENV PATH="/opt/gromacs/bin:$PATH"

# 8) Expose port for Streamlit or other web services
EXPOSE 8501

CMD ["bash", "-l"]