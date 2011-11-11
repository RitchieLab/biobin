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
   #../../data/knowledge.bio \
   #../../data/example.cfg \
   #../../data/example.vcf
   
dist_data_DATA=$(top_srcdir)/data/knowledge.bio \
	$(top_srcdir)/data/example.cfg \
	$(top_srcdir)/data/example.vcf 
   
AM_CPPFLAGS=-I$(top_srcdir)/src $(BOOST_CPPFLAGS) $(SOCI_CPPFLAGS) $(SOCI_SQLITE_CPPFLAGS) -DDATA_DIR='"$(datadir)"'
AM_LDFLAGS=$(BOOST_LDFLAGS) $(SOCI_LDFLAGS) $(SOCI_SQLITE_LDFLAGS)

biobin_LDADD=../utility/libutility.la ../biofilter/libbiofilter.la \
	../knowledge/libknowledge.la ../liftover/libliftover.la \
	../vcftools/libvcftools.la ../ldspline/libldspline.la \
	$(BOOST_REGEX_LIB) $(SOCI_LIB) $(SOCI_SQLITE_LIB)
	

