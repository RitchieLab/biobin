# ===========================================================================
#	   http://www.gnu.org/software/autoconf-archive/ax_boost_base.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_BOOST_BASE([MINIMUM-VERSION], [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   Test for the Boost C++ libraries of a particular version (or newer)
#
#   If no path to the installed boost library is given the macro searchs
#   under /usr, /usr/local, /opt and /opt/local and evaluates the
#   $BOOST_ROOT environment variable. Further documentation is available at
#   <http://randspringer.de/boost/index.html>.
#
#   This macro calls:
#
#	 AC_SUBST(BOOST_CPPFLAGS) / AC_SUBST(BOOST_LDFLAGS)
#
#   And sets:
#
#	 HAVE_BOOST
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

AC_DEFUN([AX_BOOST_BASE],
[
AC_ARG_WITH([boost],
  [AS_HELP_STRING([--with-boost@<:@=ARG@:>@],
	[use Boost library from a standard location (ARG=yes),
	 from the specified location (ARG=<path>),
	 or disable it (ARG=no)
	 @<:@ARG=yes@:>@ ])],
	[
	if test "$withval" = "no"; then
		want_boost="no"
	elif test "$withval" = "yes"; then
		want_boost="yes"
		ac_boost_path=""
	else
		want_boost="yes"
		ac_boost_path="$withval"
	fi
	],
	[want_boost="yes"])


AC_ARG_WITH([boost-libdir],
		AS_HELP_STRING([--with-boost-libdir=LIB_DIR],
		[Force given directory for boost libraries. Note that this will override library path detection, so use this parameter only if default library detection fails and you know exactly where your boost libraries are located.]),
		[
		if test -d "$withval"
		then
				ac_boost_lib_path="$withval"
		else
				AC_MSG_ERROR(--with-boost-libdir expected directory name)
		fi
		],
		[ac_boost_lib_path=""]
)

if test "x$want_boost" = "xyes"; then
	boost_lib_version_req=ifelse([$1], ,1.20.0,$1)
	boost_lib_version_req_shorten=`expr $boost_lib_version_req : '\([[0-9]]*\.[[0-9]]*\)'`
	boost_lib_version_req_major=`expr $boost_lib_version_req : '\([[0-9]]*\)'`
	boost_lib_version_req_minor=`expr $boost_lib_version_req : '[[0-9]]*\.\([[0-9]]*\)'`
	boost_lib_version_req_sub_minor=`expr $boost_lib_version_req : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
	if test "x$boost_lib_version_req_sub_minor" = "x" ; then
		boost_lib_version_req_sub_minor="0"
		fi
	WANT_BOOST_VERSION=`expr $boost_lib_version_req_major \* 100000 \+  $boost_lib_version_req_minor \* 100 \+ $boost_lib_version_req_sub_minor`
	AC_MSG_CHECKING(for boostlib >= $boost_lib_version_req)
	succeeded=no

	dnl On 64-bit systems check for system libraries in both lib64 and lib.
	dnl The former is specified by FHS, but e.g. Debian does not adhere to
	dnl this (as it rises problems for generic multi-arch support).
	dnl The last entry in the list is chosen by default when no libraries
	dnl are found, e.g. when only header-only libraries are installed!
	libsubdirs="lib"
	ax_arch=`uname -m`
	if test $ax_arch = x86_64 -o $ax_arch = ppc64 -o $ax_arch = s390x -o $ax_arch = sparc64; then
		libsubdirs="lib lib64"
	fi

	dnl first we check the system location for boost libraries
	dnl this location ist chosen if boost libraries are installed with the --layout=system option
	dnl or if you install boost with RPM
	
	std_paths="/usr /usr/local /opt /opt/local"

	if test "$ac_boost_path" != ""; then
		BOOST_CPPFLAGS="-I$ac_boost_path/include"
		CPPFLAGS_SAVED="$CPPFLAGS"
		CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
		export CPPFLAGS

		AC_REQUIRE([AC_PROG_CXX])
		AC_LANG_PUSH(C++)
		AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
			@%:@include <boost/version.hpp>
			]], [[
			#if BOOST_VERSION >= $WANT_BOOST_VERSION
			// Everything is okay
			#else
			#  error Boost version is too old
			#endif
			]])],[
			succeeded=yes
			found_system=yes
			],[
		])
		CPPFLAGS="$CPPFLAGS_SAVED"
		AC_LANG_POP(C++)

	elif test "$cross_compiling" != yes; then	
		_version=0
	
		PATH_LIST="`echo $CPLUS_INCLUDE_PATH:$C_INCLUDE_PATH | sed 's/:/ /g' | sed 's/\/include\/* */ /g'` $std_paths"
		for ac_boost_path_tmp in $PATH_LIST; do
			if test -r "$ac_boost_path_tmp/include/boost/version.hpp"; then
				BOOST_CPPFLAGS="-I$ac_boost_path_tmp/include"

				CPPFLAGS_SAVED="$CPPFLAGS"
				CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
				export CPPFLAGS

				AC_REQUIRE([AC_PROG_CXX])
				AC_LANG_PUSH(C++)
				AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
					@%:@include <boost/version.hpp>
					]], [[
					#if BOOST_VERSION >= $WANT_BOOST_VERSION
					// Everything is okay
					#else
					#  error Boost version is too old
					#endif
					]])],[
					succeeded=yes
					found_system=yes
					],[
				])
				CPPFLAGS="$CPPFLAGS_SAVED"
				AC_LANG_POP(C++)
				
				if test "x$succeeded" = "xyes"; then 				
					break
				fi
			fi
		done
		
	fi
	
	if test "x$succeeded" = "xyes"; then 
	
		BOOST_INC_PATH=`echo "$BOOST_CPPFLAGS" | sed 's/^-I//g'`
				
		# Now, get the version!
		V=`cat "$ac_boost_path_tmp/include/boost/version.hpp" | grep BOOST_VERSION | grep '#define' | grep -v HPP | sed 's/.*BOOST_VERSION//g'`
		MAJ_V=`expr $V / 100000`
		MIN_V=`expr $V / 100 % 1000`
		PATCH_V=`expr $V % 100`
	
		# Check for libraries that match the version
		# first, check relative to the include path
		if test "$ac_boost_lib_path" != ""; then
			BOOST_LDFLAGS="-L$ac_boost_lib_path"
		else
		
			found_lib="no"
					
			for libsubdir in $libsubdirs; do
				if ls "$ac_boost_path_tmp/$libsubdir/libboost_"*$MAJ_V.$MIN_V.$PATCH_V >/dev/null 2>&1; then
					BOOST_LDFLAGS="-L$ac_boost_path_tmp/$libsubdir"
					found_lib="yes"
					break
				fi
			done
			
			# Now search the library path for a matching version
			if test "x$found_lib" != "xyes"; then
				LIB_PATH_LIST="`echo $LIBRARY_PATH:$LPATH | sed 's/:/ /g' | sed 's/\/lib\(64\)*\/* */ /g'` $std_paths"
				for lib_path in $LIB_PATH_LIST; do
					for libsubdir in $libsubdirs; do
						if ls "$lib_path/$libsubdir/libboost_"*$MAJ_V.$MIN_V.$PATCH_V >/dev/null 2>&1; then
							BOOST_LDFLAGS="-L$lib_path/$libsubdir"
							found_lib="yes"
							break
						fi
					done
					if test "x$found_lib" = "xyes"; then break; fi
				done
			fi
		fi
		
		break
	fi


	#AC_LANG_POP([C++])


	if test "$succeeded" != "yes" ; then
		AC_MSG_RESULT(no)
		# execute ACTION-IF-NOT-FOUND (if present):
		ifelse([$3], , :, [$3])
	else
		AC_MSG_RESULT(yes)
		AC_SUBST(BOOST_CPPFLAGS)
		AC_SUBST(BOOST_LDFLAGS)
		if test "$BOOST_LDFLAGS" = ""; then
			AC_MSG_WARN([Could not find compiled Boost libraries, use --with-boost-libdir to force detection])
		fi
		AC_DEFINE(HAVE_BOOST,,[define if the Boost library is available])
		# execute ACTION-IF-FOUND (if present):
		ifelse([$2], , :, [$2])
	fi

	CPPFLAGS="$CPPFLAGS_SAVED"
	LDFLAGS="$LDFLAGS_SAVED"
fi

])
