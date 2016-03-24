#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(biobin, 2.2.0, software@ritchielab.psu.edu)

# Checks for --enable-debug added
AX_CHECK_ENABLE_DEBUG

AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX

RL_CONFIGURE([1.42.0],[HAVE_SQLITE HAVE_BOOST HAVE_BOOST_FILESYSTEM \
                       HAVE_BOOST_PROGRAM_OPTIONS HAVE_BOOST_SYSTEM \
                       HAVE_BOOST_IOSTREAMS HAVE_BOOST_THREAD HAVE_GSL],
	[2.7],[apsw suds])

AM_INIT_AUTOMAKE([1.9 foreign])

AC_CONFIG_SRCDIR([src/biobin/main.cpp])
AC_CONFIG_HEADER([config.h])

AC_CONFIG_FILES([Makefile
                 src/biobin/Makefile
                 src/knowledge/Makefile
                 updater/Makefile
                 updater/loki/loadPops/Makefile
                 updater/loki/loadPops/ldspline/Makefile
                 updater/loki/loadPops/ldspline/src/ldspline/Makefile
                 updater/loki/loadPops/ldspline/src/utility/Makefile])
AC_OUTPUT
