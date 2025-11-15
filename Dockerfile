FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y build-essential cmake git pkg-config libfftw3-dev libhdf5-dev libeigen3-dev libboost-all-dev libexpat1-dev libomp-dev python3 python3-pip ca-certificates && rm -rf /var/lib/apt/lists/*
WORKDIR /opt
RUN git clone -b v2025.3 https://gitlab.com/gromacs/gromacs.git && cmake -S gromacs -B gromacs-build -DGMX_BUILD_OWN_FFTW=OFF -DGMX_OPENMP=ON -DGMX_SIMD=AUTO -DCMAKE_INSTALL_PREFIX=/opt/gromacs && cmake --build gromacs-build -j && cmake --install gromacs-build
RUN git clone -b master https://github.com/votca/votca.git && cmake -S votca -B votca-build -DENABLE_CSG=ON -DBUILD_XTP=OFF -DINSTALL_CSGAPPS=ON -DWITH_GMX=ON -DCMAKE_PREFIX_PATH=/opt/gromacs -DCMAKE_INSTALL_PREFIX=/opt/votca && cmake --build votca-build -j && cmake --install votca-build
ENV PATH=/opt/votca/bin:/opt/gromacs/bin:$PATH
ENV VOTCASHARE=/opt/votca/share/votca
ENTRYPOINT ["/bin/bash"]