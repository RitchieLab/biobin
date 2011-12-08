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
   configuration.h \
   configuration.cpp \
   dataimporter.h \
   dataimporter.cpp \
   main.h \
   main.cpp \
   task.h \
   task.cpp \
   taskfilegeneration.h \
   taskfilegeneration.cpp \
   PopulationManager.h \
   PopulationManager.cpp
   #../../data/knowledge.bio \
   #../../data/example.cfg \
   #../../data/example.vcf
   
dist_data_DATA=$(top_srcdir)/data/knowledge.bio \
	$(top_srcdir)/data/example.cfg \
	$(top_srcdir)/data/example.vcf 
   
AM_CPPFLAGS=-I$(top_srcdir)/src $(BOOST_CPPFLAGS) $(SQLITE_CFLAGS) -DDATA_DIR='"$(datadir)"'
AM_LDFLAGS=$(BOOST_LDFLAGS) 

biobin_LDFLAGS=

biobin_LDADD=../utility/libutility.la \
	../knowledge/libknowledge.la \
	../vcftools/libvcftools.la \
	$(SQLITE_LIBS) \
	$(BOOST_REGEX_LIB) 

