bin_PROGRAMS = biobin

biobin_SOURCES= \
   application.h \
   application.cpp \
   Bin.h \
   Bin.cpp \
   binapplication.h \
   binapplication.cpp \
   binmanager.h \
   binmanager.cpp \
   main.h \
   main.cpp \
   task.h \
   task.cpp \
   taskfilegeneration.h \
   taskfilegeneration.cpp \
   PopulationManager.h \
   PopulationManager.cpp \
   Configuration.h \
   Configuration.cpp

EXTRA_DIST=\
	$(top_srcdir)/doc/manual.doc \
	$(top_srcdir)/updater/loadPops/buildPopulations.py
	
bin_SCRIPTS= \
	$(top_srcdir)/updater/loadPops/buildPopulations.py
   
AM_CPPFLAGS=-I$(top_srcdir)/src $(BOOST_CPPFLAGS) $(SQLITE_CFLAGS) -DDATA_DIR='"$(datadir)"'
AM_LDFLAGS=$(BOOST_LDFLAGS) 

biobin_LDADD=\
	../knowledge/libknowledge.la \
	../vcftools/libvcftools.la \
	$(SQLITE_LIBS) \
	$(BOOST_REGEX_LIB) \
	$(BOOST_FILESYSTEM_LIB) \
	$(BOOST_SYSTEM_LIB) \
	$(BOOST_PROGRAM_OPTIONS_LIB) 

