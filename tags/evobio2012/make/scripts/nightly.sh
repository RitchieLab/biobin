#! /bin/sh

# This script checks out and builds any project in the Ritchie Lab (using RL build scripts)

if test -z $PROJECT; then 
	echo "Error: You must supply the PROJECT environment variable to run this script"
	exit;
fi

RCC_GROUP_DIR=/gpfs/group/mdr23
BX_GROUP_DIR=/afs/bx.psu.edu/depot/data/ritchie_lab

if test -d $RCC_GROUP_DIR; then
	GROUP_DIR=$RCC_GROUP_DIR
elif test -d $BX_GROUP_DIR; then
	GROUP_DIR=$BX_GROUP_DIR
else
	echo "Unknown build location.. exiting"
	exit;
fi

. $GROUP_DIR/envvars.include

SVN_REPO="--username=rlb5494 --password=Build4RL --non-interactive --trust-server-cert https://ritchielab.psu.edu/svn"

DEST_DIR=$GROUP_DIR/builds/nightly/$PROJECT
INST_DIR=$GROUP_DIR/software/$PROJECT/nightly

rm -rf $DEST_DIR
svn co $SVN_REPO/projects/$PROJECT/trunk $DEST_DIR
cd $DEST_DIR
autoreconf -i
mkdir build
cd build
../configure --prefix=$INST_DIR --bindir=$INST_DIR --datadir=$INST_DIR
make -j10
make install
make -j10 distcheck
