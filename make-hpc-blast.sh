#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Incorrect number of arguments given. Need a ROOT build directory."
    exit
fi

export ROOTDIR=$1

if [ ! -d ncbi-blast-2.7.1+-src ]; then
    if [ ! -f ncbi-blast-2.7.1+-src.tar.gz ]; then
        echo "Missing the ncbi-blast-2.7.1+-src.tar.gz file in this directory."
        exit
    fi
    tar xvzf ncbi-blast-2.7.1+-src.tar.gz
fi

cp api/blast_options_handle.cpp ncbi-blast-2.7.1+-src/c++/src/algo/blast/api/
cp api/local_blast.cpp ncbi-blast-2.7.1+-src/c++/src/algo/blast/api/
cp api/traceback_stage.cpp ncbi-blast-2.7.1+-src/c++/src/algo/blast/api/

cp include/local_blast.hpp ncbi-blast-2.7.1+-src/c++/include/algo/blast/api/
cp include/traceback_stage.hpp ncbi-blast-2.7.1+-src/c++/include/algo/blast/api/

cp blast/hpc_blastp_app.cpp ncbi-blast-2.7.1+-src/c++/src/app/blast/hpc_blastp_app.cpp
cp blast/Makefile.hpc_blastp.app ncbi-blast-2.7.1+-src/c++/src/app/blast/

cp blast/hpc_blastn_app.cpp ncbi-blast-2.7.1+-src/c++/src/app/blast/
cp blast/Makefile.hpc_blastn.app ncbi-blast-2.7.1+-src/c++/src/app/blast/

cp blast/Makefile.in ncbi-blast-2.7.1+-src/c++/src/app/blast/

export CC=mpiicc
#export CC=icc
#export CC=gcc
export CXX=mpiicpc
#export CXX=icpc
#export CXX=g++
export LD=xild
export AR="xiar crs"

cd ncbi-blast-2.7.1+-src/c++/

./configure --without-3psw --with-bin-release --with-static-exe --with-mt --without-debug --without-boost --without-strip --with-build-root=$ROOTDIR

mkdir $ROOTDIR/build
cd $ROOTDIR/build
