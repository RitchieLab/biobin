SUBDIRS = ldspline

bin_PROGRAMS = pop_loader

pop_loader_SOURCES = ldsplineimporter.cpp \
	ldsplineimporter.h \
	main.cpp
	
pop_loader_LDADD = ./ldspline/src/ldspline/libldspline.la \
	$(SQLITE_LIBS) 	$(BOOST_FILESYSTEM_LIB)

AM_LDFLAGS= $(BOOST_LDFLAGS) $(BOOST_FILESYSTEM_LDFLAGS) 

AM_CPPFLAGS = -I$(srcdir)/ldspline/src -DDATA_DIR='"$(datadir)"' \
	$(BOOST_CPPFLAGS) $(BOOST_FILESYSTEM_CPPFLAGS) $(SQLITE_CFLAGS)

bin_SCRIPTS = buildPopulations.py

EXTRA_DIST = buildPopulations.py