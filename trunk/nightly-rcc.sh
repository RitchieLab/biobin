#! /bin/sh

# This script checks out and builds the biobin application

GROUP_DIR=/gpfs/group/mdr23
SVN_REPO="--username=rlb5494 --password=Build4RL --non-interactive --trust-server-cert https://ritchielab.psu.edu/svn"

DEST_DIR=$HOME/nightly/biobin
INST_DIR=$GROUP_DIR/software/biobin/nightly

rm -rf $DEST_DIR
svn co $SVN_REPO/projects/biobin/trunk $DEST_DIR
cd $DEST_DIR
autoreconf -i
mkdir build
cd build
../configure --prefix=$INST_DIR --bindir=$INST_DIR --datadir=$INST_DIR
make -j10
make install
