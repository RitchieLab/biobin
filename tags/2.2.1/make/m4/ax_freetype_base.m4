AC_DEFUN([AX_FREETYPE_BASE],
[
AC_ARG_WITH([freetype],
  [AS_HELP_STRING([--with-freetype@<:@=ARG@:>@],
    [use Boost library from a standard location (ARG=yes),
     from the specified location (ARG=<path>),
     or disable it (ARG=no)
     @<:@ARG=yes@:>@ ])],
    [
    if test "$withval" = "no"; then
        want_freetype="no"
    elif test "$withval" = "yes"; then
        want_freetype="yes"
        ac_freetype_path=""
    else
        want_freetype="yes"
        ac_freetype_path="$withval"
    fi
    ],
    [want_freetype="yes"])


AC_ARG_WITH([freetype-libdir],
        AS_HELP_STRING([--with-freetype-libdir=LIB_DIR],
        [Force given directory for freetype libraries. Note that this will override library path detection, so use this parameter only if default library detection fails and you know exactly where your freetype libraries are located.]),
        [
        if test -d "$withval"
        then
                ac_freetype_lib_path="$withval"
        else
                AC_MSG_ERROR(--with-freetype-libdir expected directory name)
        fi
        ],
        [ac_freetype_lib_path=""]
)

if test "x$want_freetype" = "xyes"; then

		freetype_compile=no
		freetype_link=no

    freetype_lib_version_req=ifelse([$1], ,2.0.9,$1)
    freetype_lib_version_req_shorten=`expr $freetype_lib_version_req : '\([[0-9]]*\.[[0-9]]*\)'`
    freetype_lib_version_req_major=`expr $freetype_lib_version_req : '\([[0-9]]*\)'`
    freetype_lib_version_req_minor=`expr $freetype_lib_version_req : '[[0-9]]*\.\([[0-9]]*\)'`
    freetype_lib_version_req_sub_minor=`expr $freetype_lib_version_req : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
    if test "x$freetype_lib_version_req_sub_minor" = "x" ; then
        freetype_lib_version_req_sub_minor="0"
        fi
    WANT_FREETYPE_VERSION=`expr $freetype_lib_version_req_major \* 100000 \+  $freetype_lib_version_req_minor \* 100 \+ $freetype_lib_version_req_sub_minor`
    
    dnl First, check for the location of the freetype libraries
    
    if test "$ac_freetype_path" != ""; then
				if test -d "$ac_freetype_path/include/freetype2" && test -r "$ac_freetype_path/include/freetype2"; then
		    		FREETYPE_CPPFLAGS="-I$ac_freetype_path/include/freetype2"
		    else
		    		FREETYPE_CPPFLAGS="-I$ac_freetype_path"
		    fi
		else
				FREETYPE_CPPFLAGS=`freetype-config --cflags`
		fi
		
		dnl Now, check the version number to see if it is correct!

    AC_MSG_CHECKING([for Freetype >= $freetype_lib_version_req])
    
    old_CPPFLAGS="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS $FREETYPE_CPPFLAGS"
    export CPPFLAGS
    AC_TRY_CPP([#include <ft2build.h> 
                #include FT_FREETYPE_H
								#if FREETYPE_MAJOR*100000+FREETYPE_MINOR*100+FREETYPE_PATCH < $WANT_FREETYPE_VERSION
								#error Freetype version too low.
								#endif
							],
        			[freetype_compile=yes],
        			[]
    )
    
    dnl Great, now check for linking!

 		FREETYPE_LIBSTR=`freetype-config --libs`   
    if test "$ac_freetype_lib_path" != ""; then
    		FREETYPE_LDFLAGS="-L$ac_freetype_lib_path"
    else
    		FREETYPE_LDFLAGS=`echo $FREETYPE_LIBSTR | sed 's/ *-l[a-zA-Z0-9]*$//'`
		fi
 		FREETYPE_LIB=`echo $FREETYPE_LIBSTR | grep -Go "\-l[a-zA-Z0-9]*$"`    
    
  	old_LDFLAGS="$LDFLAGS"
  	LDFLAGS="$LDFLAGS $FREETYPE_LDFLAGS $FREETYPE_LIB"
  	export LDFLAGS
    
    AC_LANG_PUSH([C++])    
    AC_LINK_IFELSE(
        [AC_LANG_PROGRAM([#include <ft2build.h> 
        								  #include FT_FREETYPE_H
        								  ],
            [FT_Library lib;
             FT_Init_FreeType(&lib);])],
          [freetype_link=yes],
          []
    )
    AC_LANG_POP([C++])
    
    if test "x$freetype_compile" = "xyes" -a "x$freetype_link" = "xyes"; then
        AC_MSG_RESULT(yes)
        AC_SUBST(FREETYPE_CPPFLAGS)
        AC_SUBST(FREETYPE_LDFLAGS)
        AC_SUBST(FREETYPE_LIB)
        AC_DEFINE(HAVE_FREETYPE,,[define if the freetype library is available])
        ifelse([$2], , :, [$2])
    else
    		AC_MSG_RESULT(no)
        ifelse([$3], , :, [$3])
    fi
    
    CPPFLAGS="$old_CPPFLAGS"		
    LDFLAGS="$old_LDFLAGS"
fi

])
