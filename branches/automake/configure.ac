#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(biobin, 1.0, bugs@ritchielab.psu.edu)

AM_INIT_AUTOMAKE([1.9 foreign])

AC_CONFIG_SRCDIR([src/biobin/main.cpp])
AC_CONFIG_HEADER([config.h])

m4_include([make/configure.ac])

AC_CONFIG_FILES([Makefile
                 src/biobin/Makefile
                 src/biofilter/Makefile
                 src/knowledge/Makefile
                 src/ldspline/Makefile
                 src/liftover/Makefile
                 src/utility/Makefile
                 src/vcftools/Makefile])
AC_OUTPUT
