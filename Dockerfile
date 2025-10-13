########################################################################################################################
# imsrg build stage
########################################################################################################################

FROM ubuntu:24.04 as build

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
  build-essential \
  gcc \
  libomp-dev \
  libopenblas-dev \
  liblapack-dev \
  libgsl-dev \
  libhdf5-dev \
  zlib1g-dev \
  libboost-all-dev \
  python3 \
  python3-numpy \
  cmake \
  && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN cmake . 
RUN make

FROM ubuntu:24.04 as release

RUN apt-get update && apt-get install -y \
  libomp-dev \
  libgomp1 \
  libopenblas-dev \
  libgsl-dev \
  python3 \
  python3-numpy \
  python3-sympy \
  python3-pandas \
  && apt-get clean && rm -rf /var/lib/apt/lists/*

COPY pyIMSRG.cpython-312-x86_64-linux-gnu.so /usr/local/lib/

ENV PYTHONPATH="$PYTHONPATH:/usr/local/lib/"
ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"
  