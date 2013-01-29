noinst_LTLIBRARIES = libldspline.la

libldspline_la_LIBADD = \
	../utility/libld_utility.la

AM_CPPFLAGS = \
	-I$(srcdir)/.. 

libldspline_la_SOURCES = \
	ldspline.h \
	locuslookup.h \
	snpspline.h \
	timestamp.h \
	ldspline.cpp \
	locuslookup.cpp \
	snpspline.cpp
	
bin_PROGRAMS = ldspline

ldspline_SOURCES = 	main.cpp
	
ldspline_LDADD = ../utility/libld_utility.la libldspline.la
	
