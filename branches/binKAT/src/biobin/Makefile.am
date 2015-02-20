bin_PROGRAMS = biobin

biobin_SOURCES= \
   Bin.h \
   Bin.cpp \
   binapplication.h \
   binapplication.cpp \
   binmanager.h \
   binmanager.cpp \
   main.h \
   main.cpp \
   PopulationManager.h \
   PopulationManager.cpp \
   Configuration.h \
   Configuration.cpp \
   util/ICompressedFile.h \
   util/ICompressedFile.cpp \
   util/string_ref.hpp \
   util/string_ref_fwd.hpp \
   util/Phenotype.h \
   tests/Test.h \
   tests/Test.cpp \
   tests/TestFactory.h \
   tests/TestFactory.cpp \
   tests/LinearRegression.h \
   tests/LinearRegression.cpp \
   tests/MatrixUtils.h \
   tests/MatrixUtils.cpp \
   tests/SKATLinear.h \
   tests/SKATLinear.cpp \
   tests/qfc.h \
   tests/qfc.cpp \
   tests/Wilcoxon.h \
   tests/Wilcoxon.cpp

EXTRA_DIST=$(top_srcdir)/doc/manual.odt 
 
AM_CPPFLAGS=-I$(top_srcdir)/src $(BOOST_CPPFLAGS) $(SQLITE_CFLAGS) $(GSL_CFLAGS) -DDATA_DIR='"$(datadir)"'
AM_LDFLAGS=$(BOOST_LDFLAGS) $(GSL_LDFLAGS)


biobin_LDADD=\
	../knowledge/libknowledge.la \
	$(SQLITE_LIBS) \
	$(BOOST_REGEX_LIB) \
	$(BOOST_FILESYSTEM_LIB) \
	$(BOOST_SYSTEM_LIB) \
	$(BOOST_PROGRAM_OPTIONS_LIB) \
	$(BOOST_IOSTREAMS_LIB) \
	$(BOOST_THREAD_LIB) \
	$(GSL_LIBS)

