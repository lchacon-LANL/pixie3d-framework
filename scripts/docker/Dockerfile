# Use latest Ubuntu as base image
FROM ubuntu:bionic

# prevent user interaction
ARG DEBIAN_FRONTEND=noninteractive

# Install any needed packages
RUN apt-get update \
    && apt-get install -y apt-transport-https ca-certificates gnupg software-properties-common wget \
    && wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null \
    && apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main' \
    && apt-get update \
	&& apt-get install -y \
		gcc \
		g++ \
		gfortran \
		cmake \
		git \
		python \
		python3 \
		python3-pip \
		openmpi-bin \
		libopenmpi-dev \
		libblas-dev \
		liblapack-dev \
		autoconf \
	&& rm -rf /var/lib/apt/lists/*

# Define environment variable
ENV CC /usr/bin/gcc
ENV CXX /usr/bin/g++
ENV FC /usr/bin/gfortran
ENV F77 /usr/bin/gfortran

ENV OMPI_CC /usr/bin/gcc
ENV OMPI_CXX /usr/bin/g++
ENV OMPI_FC /usr/bin/gfortran

ENV MPICC /usr/bin/mpicc
ENV MPICXX /usr/bin/mpicxx
ENV MPIFC /usr/bin/mpifort

# install fortran format
RUN pip3 install --upgrade fprettify cmake-format

# using tmp dir
WORKDIR /tmp

# set install folder
ENV LIB_INSTALL_DIR /opt/


# build petsc
COPY petsc-3.0.0-p7.tgz .
RUN tar -xaf petsc-3.0.0-p7.tgz \
	&& cd petsc-3.0.0-p7 \
	&& ./configure --prefix=$LIB_INSTALL_DIR --with-cc=mpicc --with-fc=mpifort --with-shared -- with-fortran --with-debugging=1 \
	&& export PETSC_ARCH=linux-gnu-c-debug \
	&& export PETSC_DIR=$(pwd) \
	&& make -j \
	&& make install \
	&& cd .. \
	&& rm -r petsc-3.0.0-p7*

ENV PETSC_DIR /opt/

# build adios
COPY adios-1.10.0.tar.gz .
RUN tar -xaf adios-1.10.0.tar.gz \
	&& cd adios-1.10.0 \
	&& mkdir build \
	&& cd build \
	&& ../configure --prefix=$LIB_INSTALL_DIR --enable-dependency-tracking --without-netcdf --without-nc4par --without-hdf5 \
	&& make \
	&& make install \
	&& cd ../../ \
	&& rm -r adios-1*

# build sprng
COPY sprng5.tar.gz .
RUN tar -xaf sprng5.tar.gz \
	&& cd sprng5 \
	&& ./configure --prefix=$LIB_INSTALL_DIR --with-fortran \ 
	&& make \
	&& make install \
	&& cp include/* $LIB_INSTALL_DIR/include/ \
	&& cd ../ \
	&& rm -r sprng5*

