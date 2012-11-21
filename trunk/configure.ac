#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(biobin, 2.0.3, software@ritchielab.psu.edu)

AM_INIT_AUTOMAKE([1.9 foreign])

AC_CONFIG_SRCDIR([src/biobin/main.cpp])
AC_CONFIG_HEADER([config.h])

RL_CONFIGURE([1.46.0],[HAVE_BOOST HAVE_BOOST_THREAD HAVE_BOOST_FILESYSTEM HAVE_SQLITE HAVE_BOOST_PROGRAM_OPTIONS HAVE_BOOST_SYSTEM],
	[2.7],[apsw suds])

AC_CONFIG_FILES([Makefile
                 src/biobin/Makefile
                 src/knowledge/Makefile
                 src/vcftools/Makefile
                 loki/Makefile
                 loki/loadPops/Makefile
                 loki/loadPops/ldspline/Makefile
                 loki/loadPops/ldspline/src/ldspline/Makefile
                 loki/loadPops/ldspline/src/utility/Makefile])
AC_OUTPUT
