noinst_LTLIBRARIES = libldspline.la

libldspline_la_LIBADD = \
	../utility/libld_utility.la $(BOOST_IOSTREAMS_LIB)

AM_CPPFLAGS = \
	-I$(srcdir)/..
	
AM_LDFLAGS=$(BOOST_LDFLAGS)

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

ldspline_LDADD = ./libldspline.la ../utility/libld_utility.la $(BOOST_IOSTREAMS_LIB)
	
	
