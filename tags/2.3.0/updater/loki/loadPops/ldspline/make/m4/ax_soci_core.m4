# ===========================================================================
#	   http://www.gnu.org/software/autoconf-archive/ax_soci_base.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_SOCI_CORE([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   Test for the Boost C++ libraries of a particular version (or newer)
#
#   If no path to the installed soci library is given the macro searchs
#   under /usr, /usr/local, /opt and /opt/local and evaluates the
#   $SOCI_ROOT environment variable. Further documentation is available at
#   <http://randspringer.de/soci/index.html>.
#
#   This macro calls:
#
#	 AC_SUBST(SOCI_CPPFLAGS) / AC_SUBST(SOCI_LDFLAGS)
#
#   And sets:
#
#	 HAVE_SOCI
#
# LICENSE
#
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2009 Peter Adolphs
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 20

AC_DEFUN([AX_SOCI_CORE],
[
AC_ARG_WITH([soci],
  [AS_HELP_STRING([--with-soci@<:@=ARG@:>@],
	[use Boost library from a standard location (ARG=yes),
	 from the specified location (ARG=<path>),
	 or disable it (ARG=no)
	 @<:@ARG=yes@:>@ ])],
	[
	if test "$withval" = "no"; then
		want_soci="no"
	elif test "$withval" = "yes"; then
		want_soci="yes"
		ac_soci_path=""
	else
		want_soci="yes"
		ac_soci_path="$withval"
	fi
	],
	[want_soci="yes"])


AC_ARG_WITH([soci-libdir],
		AS_HELP_STRING([--with-soci-libdir=LIB_DIR],
		[Force given directory for soci libraries. Note that this will override library path detection, so use this parameter only if default library detection fails and you know exactly where your soci libraries are located.]),
		[
		if test -d "$withval"
		then
			ac_soci_lib_path="$withval"
		else
			AC_MSG_ERROR(--with-soci-libdir expected directory name)
		fi
		],
		[ac_soci_lib_path=""]
)


if test "x$want_soci" = "xyes"; then
	AC_CHECK_LIB([dl], [main])
	AC_MSG_CHECKING(for SOCI core libraries)
	AC_REQUIRE([AC_PROG_CXX])
	AC_LANG_PUSH(C++)
	succeeded=no

	dnl Check for lib64 in appropriate systems
	
	libsubdirs="lib"
	ax_arch=`uname -m`
	if test $ax_arch = x86_64 -o $ax_arch = ppc64 -o $ax_arch = s390x -o $ax_arch = sparc64; then
		libsubdirs="lib64 lib"
	fi

	dnl List all the standard search paths here
	std_paths="/usr/local /usr /opt /opt/local"
	
	if test "$ac_soci_path" != ""; then
		if test -d "$ac_soci_path/include/soci" && test -r "$ac_soci_path/include/soci"; then
			SOCI_CPPFLAGS="-I$ac_soci_path/include/soci"
		else
			SOCI_CPPFLAGS="-I$ac_soci_path/include"
		fi
		
		CPPFLAGS_SAVED="$CPPFLAGS"
		CPPFLAGS="$CPPFLAGS $SOCI_CPPFLAGS"
		export CPPFLAGS  
		AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[@%:@include <soci.h>]], [[]])],[include_succeed=yes],[])
		
	else
		dnl First, check to see if the header is right there in the path!
	
		AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[@%:@include <soci.h>]], [[]])],[include_succeed=yes],[])

		if test "x$include_succeed" != "xyes"; then
			dnl Check the CPLUS_INCLUDE_PATH as well as /usr /usr/local /opt /opt/local				
			PATH_LIST="`echo $CPLUS_INCLUDE_PATH | sed 's/:/ /g' | sed 's/\/include\/* */ /g'` $std_paths"
			for ac_soci_path_tmp in $PATH_LIST; do
				if test -r "$ac_soci_path_tmp/include/soci/soci.h"; then
					SOCI_CPPFLAGS="-I$ac_soci_path_tmp/include/soci"
					break;
				elif test -r "$ac_soci_path_tmp/include/soci.h"; then
					SOCI_CPPFLAGS="-I$ac_soci_path_tmp/include"
					break;
				fi
			done
						  
			CPPFLAGS_SAVED="$CPPFLAGS"
			CPPFLAGS="$CPPFLAGS $SOCI_CPPFLAGS"
			export CPPFLAGS  
			AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[@%:@include <soci.h>]], [[]])],[include_succeed=yes],[])
			
		fi
	fi



	dnl overwrite ld flags if we have required special directory with
	dnl --with-soci-libdir parameter
	SOCI_LDFLAGS=""
	SOCI_LIB="-lsoci_core"
	
	if test "$ac_soci_lib_path" != ""; then
		SOCI_LDFLAGS="-L$ac_soci_lib_path"
			
		LDFLAGS_SAVED="$LDFLAGS"
		LDFLAGS="$LDFLAGS $SOCI_LDFLAGS $SOCI_LIB"
		export LDFLAGS
	

		dnl Check if we can link against the SOCI libraries
		AC_LINK_IFELSE(
			[AC_LANG_PROGRAM([#include <soci.h>],
				[soci::session dummy])],
			  [link_succeed=yes],
			  []
		)
		
	else
		if test "$ac_soci_path" != ""; then
			for libsubdir in $libsubdirs ; do
				if ls "$ac_soci_path/$libsubdir/libsoci_"* >/dev/null 2>&1 ; then 
					SOCI_LDFLAGS="-L$ac_soci_path/$libsubdir"
					break;
				fi
			done
			
			LDFLAGS_SAVED="$LDFLAGS"
			LDFLAGS="$LDFLAGS $SOCI_LDFLAGS $SOCI_LIB"
			export LDFLAGS
	

			dnl Check if we can link against the SOCI libraries
			AC_LINK_IFELSE(
				[AC_LANG_PROGRAM([#include <soci.h>],
					[soci::session dummy])],
				  [link_succeed=yes],
				  []
			)
		fi
		
		if test "x$link_succeed" != "xyes"; then

			LDFLAGS_SAVED="$LDFLAGS"
			LDFLAGS="$LDFLAGS $SOCI_LIB"
			export LDFLAGS
	

			dnl Check if we can link against the SOCI libraries
			AC_LINK_IFELSE(
				[AC_LANG_PROGRAM([#include <soci.h>],
					[soci::session dummy])],
				  [link_succeed=yes],
				  []
			)
	  
			if test "x$link_succeed" != "xyes"; then
			
				dnl Check the LPATH as well as the standards above
				PATH_LIST="`echo $LPATH | sed 's/:/ /g' | sed 's/\/lib\(64\)*\/* */ /g'` $std_paths"
				for ac_soci_path_tmp in $PATH_LIST ; do
					for libsubdir in $libsubdirs ; do
						if ls "$ac_soci_path_tmp/$libsubdir/libsoci_"* >/dev/null 2>&1 ; then 
							SOCI_LDFLAGS="-L$ac_soci_path_tmp/$libsubdir"
							found_lib="yes"
							break;
						fi
					done
					if test "x$found_lib" = "xyes"; then break; fi
				done
				LDFLAGS="$LDFLAGS $SOCI_LDFLAGS $SOCI_LIB"
				export LDFLAGS
				dnl Check if we can link against the SOCI libraries
				AC_LINK_IFELSE(
					[AC_LANG_PROGRAM([#include <soci.h>],
						[soci::session dummy])],
		  			[link_succeed=yes],
		  			[]
				)
			fi
		fi
	fi
	

	  
	if test "x$include_succeed" = "xyes" -a "x$link_succeed" = "xyes"; then
		AC_MSG_RESULT(yes)
		AC_SUBST(SOCI_CPPFLAGS)
		AC_SUBST(SOCI_LDFLAGS)
		AC_SUBST(SOCI_LIB)
		AC_DEFINE(HAVE_SOCI,,[define if the soci library is available])
		ifelse([$1], , :, [$1])
	else
		AC_MSG_RESULT(no)
		ifelse([$2], , :, [$2])
	fi
	
	
	AC_LANG_POP([C++])
	CPPFLAGS="$CPPFLAGS_SAVED"
	LDFLAGS="$LDFLAGS_SAVED"
fi

])
