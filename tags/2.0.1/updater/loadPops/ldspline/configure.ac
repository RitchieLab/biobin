# Ritchie Lab configure.ac template

AC_PREREQ(2.59)
AC_INIT(ldspline, 1.1.0, software@ritchielab.psu.edu)

AM_INIT_AUTOMAKE([1.8 foreign])

AC_CONFIG_SRCDIR([src/ldspline/main.cpp])
AC_CONFIG_HEADER([config.h])

RL_CONFIGURE([1.29],[HAVE_BOOST])

AC_CONFIG_FILES([Makefile
                 src/ldspline/Makefile
                 src/utility/Makefile
                ])
AC_OUTPUT
