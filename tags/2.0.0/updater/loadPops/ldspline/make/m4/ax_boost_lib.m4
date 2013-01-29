# This macro checks for the existence and linkability of a given boost
# library.  the macro takes up to 4 arguments (with the 1st 2 being mandatory)
# The arguments are:
#   1 - the name of the libarary, as it will be linked (lowercase + underscore)
#   2 - A short C++ line that will fail on compile and link if the library is not found
#   3 - (Optional) Extra libraries to include
#   4 - (Optional) A header that defines the library.  If not given, defaults to
#		boost/$1.hpp

# Note: This macro is designed to be a generalization of the 
# AX_BOOST_* libraries used previously

AC_DEFUN([AX_BOOST_LIB],
[
	
	LIB_NAME=$1
	LIB_NAME_DASH="`echo $LIB_NAME | sed 's/_/-/g'`"
	LIB_NAME_UPPER="`echo $LIB_NAME | tr [a-z] [A-Z]`"

	
	AC_ARG_WITH(m4_bpatsubst($1,[_],[-]),
		AS_HELP_STRING([--with-boost-]m4_bpatsubst($1,[_],[-])[@<:@=special-lib@:>@],
			[use the ]$1[ library from boost - it is possible to specify a certain library for the linker
						e.g. --with-boost-]m4_bpatsubst($1,[_],[-])[=]$1[-gcc-mt ]),
		[
		if test "$withval" = "no"; then
			want_lib="no"
		elif test "$withval" = "yes"; then
			want_lib="yes"
			ax_boost_lib=""
		else
			want_lib="yes"
			ax_boost_lib="$withval"
		fi
		],
		[want_lib="yes"]
	)
	
	if test "x$want_lib" = "xyes"; then
	
		HDR_INCLUDE="boost/$1.hpp"
		if test "x"$4 != "x"; then
			HDR_INCLUDE="$4"
		fi
	
		AC_MSG_CHECKING([whether the boost::]$1[ library is available])
	
		AC_REQUIRE([AX_BOOST_BASE])
		AC_LANG_PUSH([C++])
		
		CPPFLAGS_SAVE="$CPPFLAGS"
		CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
		export CPPFLAGS
		
		LDFLAGS_SAVE="$LDFLAGS"
		
		
		AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[@%:@include <$HDR_INCLUDE>]], [[$2]]),
			[lib_compile=yes],
			[lib_compile=no])
			
		if test "x$lib_compile" != "xno"; then
			# Try to find the location and name of the library
			BOOSTLIBDIR=`echo $BOOST_LDFLAGS | sed -e 's/@<:@^\/@:>@*//'`
			LIB_LIST=`ls $BOOSTLIBDIR/*boost_$1* | sed 's/.*\///g' | sed 's/\..*//g' | sed 's/^lib//g' | sort -u`
			
			LIB_FINAL=""
			
			for lib_cand in $LIB_LIST; do
				LDFLAGS="$LDFLAGS_SAVE $BOOST_LDFLAGS -l$lib_cand $3"
				export LDFLAGS
				
				AC_LINK_IFELSE(AC_LANG_PROGRAM([[@%:@include <$HDR_INCLUDE>]], [[$2]]),
					[lib_link=yes],
					[lib_link=no])
				
				if test "x$lib_link" = "xyes"; then
					LIB_FINAL="$lib_cand"
					break
				fi
			done	
			
		fi
		
		CPPFLAGS="$CPPFLAGS_SAVE"
		LDFLAGS="$LDFLAGS_SAVE"
		
		if test "x$lib_compile" = "xyes" -a "x$lib_link" = "xyes"; then
			AC_MSG_RESULT(yes)
			AC_DEFINE([HAVE_BOOST_]m4_toupper($1),,[define if the boost::]$1[ library is available])
			AC_SUBST([BOOST_]m4_toupper($1)[_LIB],"-l$LIB_FINAL $3")
		else
			AC_MSG_RESULT(no)
		fi
	fi

	
])
