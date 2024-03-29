# ===========================================================================
#         http://autoconf-archive.cryp.to/ax_check_enable_debug.html
# ===========================================================================
#
# SYNOPSIS
#
#   Check for the presence of an --enable-debug option to configure and
#   allow/avoid compiled debugging flags appropriately.
#
#   AX_CHECK_ENABLE_DEBUG([enable by default=yes/info/profile/no])
#
# DESCRIPTION
#
#   Check for the presence of an --enable-debug option to configure, with the
#   specified default value used when the option is not present.  Return the
#   value in the variable $ax_enable_debug.
#
#   Specifying 'yes' adds '-g -O0' to the compilation flags for all languages.
#   Specifying 'info' adds '-g' to the compilation flags.  Specifying 'profile'
#   adds '-g -pg' to the compilation flags and '-pg' to the linking flags.
#   Otherwise, nothing is added.  Define NDEBUG if debug disabled.  If debug
#   not enabled, ensure AC_PROG_* will not add debugging flags.  Should be
#   invoked prior to any AC_PROG_* compiler checks.
#
# LAST MODIFICATION
#
#   2011-03-09
#
# COPYLEFT
#
#   Copyright (c) 2011 Rhys Ulerich <rhys.ulerich@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_CHECK_ENABLE_DEBUG],[
    AC_BEFORE([$0],[AC_PROG_CC])dnl
    AC_BEFORE([$0],[AC_PROG_CXX])dnl
    AC_BEFORE([$0],[AC_PROG_F77])dnl
    AC_BEFORE([$0],[AC_PROG_FC])dnl

    AC_MSG_CHECKING(whether to enable debugging)
    m4_define(ax_enable_debug_default,[m4_tolower(m4_normalize(ifelse([$1],,[no],[$1])))])
    AC_ARG_ENABLE(debug,
        [AS_HELP_STRING([--enable-debug]@<:@=ax_enable_debug_default@:>@,[compile with debugging; one of yes/info/profile/no])],
        [],enable_debug=ax_enable_debug_default)
    if test "x$enable_debug" = "xyes" || test "x$enable_debug" = "x"; then
        AC_MSG_RESULT(yes)
        CFLAGS="${CFLAGS} -g -pg -Wall"
        CXXFLAGS="${CXXFLAGS} -g -O0 -Wall"
        FFLAGS="${FFLAGS} -g -O0 -Wall"
        FCFLAGS="${FCFLAGS} -g -O0 -Wall"
        OBJCFLAGS="${OBJCFLAGS} -g -O0 -Wall"
    else
        if test "x$enable_debug" = "xinfo"; then
            AC_MSG_RESULT(info)
            CFLAGS="${CFLAGS} -g"
            CXXFLAGS="${CXXFLAGS} -g"
            FFLAGS="${FFLAGS} -g"
            FCFLAGS="${FCFLAGS} -g"
            OBJCFLAGS="${OBJCFLAGS} -g"
        elif test "x$enable_debug" = "xprofile"; then
            AC_MSG_RESULT(profile)
            CFLAGS="${CFLAGS} -g -pg"
            CXXFLAGS="${CXXFLAGS} -g -pg"
            FFLAGS="${FFLAGS} -g -pg"
            FCFLAGS="${FCFLAGS} -g -pg"
            OBJCFLAGS="${OBJCFLAGS} -g -pg"
            LDFLAGS="${LDFLAGS} -pg"
        else
            AC_MSG_RESULT(no)
            dnl Ensure AC_PROG_CC/CXX/F77/FC/OBJC will not enable debug flags
            dnl by setting any unset environment flag variables
            if test "x${CFLAGS+set}" != "xset"; then
                CFLAGS="-O2"
            fi
            if test "x${CXXFLAGS+set}" != "xset"; then
                CXXFLAGS="-O2"
            fi
            if test "x${FFLAGS+set}" != "xset"; then
                FFLAGS="-O2"
            fi
            if test "x${FCFLAGS+set}" != "xset"; then
                FCFLAGS="-O2"
            fi
            if test "x${OBJCFLAGS+set}" != "xset"; then
                OBJCFLAGS="-O2"
            fi
        fi
        dnl assert.h is a NOP if NDEBUG is defined, so define it.
        AC_DEFINE(NDEBUG,,[Define if debugging is disabled])
    fi
    ax_enable_debug=$enable_debug
])
