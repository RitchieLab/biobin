bin_PROGRAMS = biobin

biobin_SOURCES=binapplication.cpp \
   binmanager.cpp \
   configuration.cpp \
   dataimporter.cpp \
   genotypestorage.cpp \
   main.cpp \
   individual.cpp \
   taskbincollapse.cpp \
   taskfilegeneration.cpp \
   binapplication.h \
   binmanager.h \
   configuration.h \
   dataimporter.h \
   genotypestorage.h \
   individual.h \
   main.h \
   taskbincollapse.h \
   taskfilegeneration.h \
   timestamp.h
   
AM_CPPFLAGS=-I$(top_srcdir)/src $(BOOST_CPPFLAGS) $(SOCI_CPPFLAGS) $(SOCI_SQLITE_CPPFLAGS)
AM_LDFLAGS=$(BOOST_LDFLAGS) $(SOCI_LDFLAGS) $(SOCI_SQLITE_LDFLAGS)

biobin_LDADD=../utility/libutility.la ../biofilter/libbiofilter.la \
	../knowledge/libknowledge.la ../liftover/libliftover.la \
	../vcftools/libvcftools.la ../ldspline/libldspline.la \
	$(BOOST_REGEX_LIB) $(SOCI_LIB) $(SOCI_SQLITE_LIB)
	

