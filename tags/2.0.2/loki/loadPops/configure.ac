# Ritchie Lab configure.ac template

AC_PREREQ(2.59)
AC_INIT(pop_loader, 1.1.0, software@ritchielab.psu.edu)

AM_INIT_AUTOMAKE([1.8 foreign])

AC_CONFIG_SRCDIR([main.cpp])
AC_CONFIG_HEADER([config.h])

RL_CONFIGURE([1.29],[HAVE_BOOST HAVE_BOOST_FILESYSTEM HAVE_SQLITE])

AC_CONFIG_FILES([Makefile
                 ldspline/Makefile
                 ldspline/src/ldspline/Makefile
                 ldspline/src/utility/Makefile
                ])
AC_OUTPUT
