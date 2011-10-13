#! /bin/sh

# This script checks out and builds the biobin application

DEST_DIR=$HOME/nightly/biobin
INST_DIR=/afs/bx.psu.edu/depot/data/ritchie_lab/software/nightly/

rm -rf $DEST_DIR
svn co svn+ssh://coltrane/afs/bx.psu.edu/depot/data/ritchie_lab/svn/projects/biobin/trunk $DEST_DIR
cd $DEST_DIR
autoreconf -i
mkdir build
cd build
../configure --prefix=$INST_DIR --program-suffix="-nightly"
make -j10
make install
