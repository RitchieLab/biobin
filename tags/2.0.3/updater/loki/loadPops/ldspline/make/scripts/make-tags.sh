#! /bin/bash -l

# This script checks out and builds any project in the Ritchie Lab (using RL build scripts)

if test -z $PROJECT; then 
	echo "Error: You must supply the PROJECT environment variable to run this script"
	exit;
fi

if [ -n "$FORCE" ] && [ "$FORCE" -eq 0 ]; then
		FORCE=""
fi

SVN_USER="--username=rlb5494 --password=Build4RL --non-interactive --trust-server-cert"
SVN_REPO="https://ritchielab.psu.edu/svn"

# Establish the group environment
# NOTE: the SVN_CO_DIR should be somewhere we have access to
SVN_CO_DIR=`mktemp -d`
svn co -q $SVN_USER $SVN_REPO/projects/infrastructure/trunk $SVN_CO_DIR --depth empty >/dev/null 2>&1
svn up $SVN_USER $SVN_CO_DIR/ritchielab.bash_profile >/dev/null 2>&1
. $SVN_CO_DIR/ritchielab.bash_profile >/dev/null 2>&1
rm -rf $SVN_CO_DIR

if test -z "$RITCHIELAB_GROUP_DIR"; then
	echo "Unable to identify Group Storage Directory"
	exit;
fi


for tag in $(svn ls $SVN_USER $SVN_REPO/projects/$PROJECT/tags); do

	tag=$(echo $tag | sed 's|/\s*$||')
	DEST_DIR=$RITCHIELAB_GROUP_DIR/builds/$PROJECT/$tag/$RITCHIELAB_OS
	INST_DIR=$RITCHIELAB_GROUP_DIR/software/$PROJECT/$tag

	rm -rf $DEST_DIR
	svn co -q $SVN_USER $SVN_REPO/projects/$PROJECT/tags/$tag $DEST_DIR
	
	# Only build if I can (uses automake) and it either doesn't exist or I'm forcing
	if test \( -n "$FORCE" -o ! -d "$INST_DIR/$RITCHIELAB_OS" \) -a -f $DEST_DIR/configure.ac; then
	
		echo "Building $PROJECT v. $tag"
	
		cd $DEST_DIR
		autoreconf -i >/dev/null 2>&1
		mkdir build
		cd build
		../configure --prefix=$INST_DIR/common --exec-prefix=$INST_DIR/$RITCHIELAB_OS --datadir=$RITCHIELAB_GROUP_DIR/datasets/loki --disable-loki >/dev/null 2>&1
		make -j10 >/dev/null 2>&1
		make install-exec >/dev/null 2>&1
		if test `make -j10 distcheck 2>&1 | tail -4 | grep -c "ready for distribution"` -lt 1; then
			make clean >/dev/null 2>&1
			echo "$PROJECT v. $tag build failed"
		else
			# At this point, the build succeeded
			for f in $(ls $INST_DIR/$RITCHIELAB_OS/bin); do
				FN=$(echo $f | sed 's/\.[^\.]*$//g')
				SUFF=$(echo $f | sed "s/^$FN//g")
				ln -fs $INST_DIR/$RITCHIELAB_OS/bin/$f $RITCHIELAB_GROUP_DIR/usr/software/$RITCHIELAB_OS/bin/$FN-$tag$SUFF
			done
		fi
	
	fi

done
