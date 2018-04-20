#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Incorrect number of arguments given. Need a ROOT build directory and build environment (HOST or MIC) and version (0.5 or 1.0)"
    exit
fi

if [ "$2" == "HOST" ]; then
    echo "Building for the Xeon host..."
elif [ "$2" == "MIC" ]; then
    if [ -z "$NCBI_DATATOOL_PATH" ]; then
        echo "The NCBI_DATATOOL_PATH variable is not set."
        exit
    fi

    echo "Building for the MIC coprocessor in native mode..."
else
    echo "Enter either HOST or MIC in the second argument to select the correct build environment"
    exit
fi

if [ "$3" == "1.0" ]; then
    echo "Building version 1.0..."
elif [ "$3" == "0.5" ]; then
    echo "Building version 0.5..."
else
    echo "Enter either 0.5 or 1.0 for the build version."
    exit 0
fi

export ROOTDIR=$1
export BUILDENV=$2
export VERSION=$3

# Append the ROOTDIR directory with the version number.
export ROOTDIR=$1-$3

if [ ! -d ncbi-blast-2.2.31+-src ]; then
    tar xvzf ncbi-blast-2.2.31+-src.tar.gz
fi

if [ "`grep optimization_level ncbi-blast-2.2.31+-src/c++/src/algo/blast/core/blast_gapalign.c`" == "" ]; then
    echo "Grep didn't see that the pragma had been added so doing it now..."
    awk -v n=3401 -v s="#pragma intel optimization_level 1" 'NR == n {print s} {print}' ncbi-blast-2.2.31+-src/c++/src/algo/blast/core/blast_gapalign.c > file.new
    mv file.new ncbi-blast-2.2.31+-src/c++/src/algo/blast/core/blast_gapalign.c
fi

if [ "`grep __MIC__ ncbi-blast-2.2.31+-src/c++/include/corelib/ncbifloat.h`" == "" ]; then
    echo "Grep didn't see the MIC flag so adding it now..."
    sed -i -e 's/defined(_GLIBCXX_CONSTEXPR)/defined(_GLIBCXX_CONSTEXPR) \&\& !defined(__MIC__)/g' ncbi-blast-2.2.31+-src/c++/include/corelib/ncbifloat.h
fi

BLASTPFILE=blast/hpc_blastp_app.cpp

if [ "$VERSION" == "1.0" ]; then
    if [ "`grep "UPDATE 1" $BLASTPFILE`" == "" ]
    then
	LINENUM=`grep -n "#define UPDATE" $BLASTPFILE | cut -f 1 -d":"`
	sed -i "${LINENUM}s/.*/#define UPDATE 1           \/\/enable thread level agglomeration of results/" $BLASTPFILE
    fi

    cp api/blast_options_handle.cpp ncbi-blast-2.2.31+-src/c++/src/algo/blast/api/
    cp api/local_blast.cpp ncbi-blast-2.2.31+-src/c++/src/algo/blast/api/
    cp api/traceback_stage.cpp ncbi-blast-2.2.31+-src/c++/src/algo/blast/api/

    cp include/local_blast.hpp ncbi-blast-2.2.31+-src/c++/include/algo/blast/api/
    cp include/traceback_stage.hpp ncbi-blast-2.2.31+-src/c++/include/algo/blast/api/

    cp $BLASTPFILE ncbi-blast-2.2.31+-src/c++/src/app/blast/hpc_blastp_app.cpp

    cp blast/Makefile.hpc_blastp.app ncbi-blast-2.2.31+-src/c++/src/app/blast/
    cp blast/Makefile.in ncbi-blast-2.2.31+-src/c++/src/app/blast/
else
    if [ "`grep "UPDATE 0" $BLASTPFILE`" == "" ]
    then
	LINENUM=`grep -n "#define UPDATE" $BLASTPFILE | cut -f 1 -d":"`
	sed -i "${LINENUM}s/.*/#define UPDATE 0           \/\/disable thread level agglomeration of results/" $BLASTPFILE
    fi

    cp api/blast_options_handle.cpp ncbi-blast-2.2.31+-src/c++/src/algo/blast/api/
    cp api/local_blast.cpp.ORIG ncbi-blast-2.2.31+-src/c++/src/algo/blast/api/local_blast.cpp
    cp api/traceback_stage.cpp.ORIG ncbi-blast-2.2.31+-src/c++/src/algo/blast/api/traceback_stage.cpp
    cp include/local_blast.hpp.ORIG ncbi-blast-2.2.31+-src/c++/include/algo/blast/api/local_blast.hpp
    cp include/traceback_stage.hpp.ORIG ncbi-blast-2.2.31+-src/c++/include/algo/blast/api/traceback_stage.hpp

    cp $BLASTPFILE ncbi-blast-2.2.31+-src/c++/src/app/blast/hpc_blastp_app.cpp

    cp blast/Makefile.hpc_blastp.app ncbi-blast-2.2.31+-src/c++/src/app/blast/

    cp blast/Makefile.in ncbi-blast-2.2.31+-src/c++/src/app/blast/
fi

export CC=mpiicc
export CXX=mpiicpc
export LD=xild
export AR="xiar crs"

cd ncbi-blast-2.2.31+-src/c++/

./configure --without-3psw --with-bin-release --with-static-exe --with-mt --without-debug --without-boost --without-strip --with-build-root=$ROOTDIR

cd $ROOTDIR/build

#Edit the Makefile.mk

# silence new Intel Compiler warning.
sed -i -e 's/openmp/qopenmp/g' Makefile.mk

if [ "$BUILDENV" == "HOST" ]; then
    sed -i -e 's/-O$/-O2 /g' Makefile.mk
    sed -i -e 's/O2/O3/g' Makefile.mk
    sed -i -e 's/-O3 -axSSSE3/-O3 -xHost/g' Makefile.mk
else
    sed -i -e 's/-O$/-O2 /g' Makefile.mk
    sed -i -e 's/O2/O3 -mmic/g' Makefile.mk
    sed -i -e 's/axSSSE3/mmic/g' Makefile.mk
fi

#make -j8 -k all_r

# Set the path to the datatool. Used for building for the MICs.
#if [ "$BUILDENV" == "HOST" ]; then
#    cd ../
#    export NCBI_DATATOOL_PATH=$PWD
#fi
