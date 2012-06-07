noinst_LTLIBRARIES=libknowledge.la

libknowledge_la_SOURCES=\
	RegionCollection.h \
	RegionCollection.cpp \
	RegionCollectionSQLite.h \
	RegionCollectionSQLite.cpp \
	Region.h  \
	Region.cpp \
	Allele.h \
	Allele.cpp \
	Locus.h \
	Locus.cpp \
	Group.h \
	Group.cpp \
	GroupCollection.h \
	GroupCollection.cpp \
	GroupCollectionSQLite.h \
	GroupCollectionSQLite.cpp \
	Information.h \
	InformationSQLite.h \
	InformationSQLite.cpp \
	Configuration.h \
	Configuration.cpp \
	liftover/BuildConversion.h \
	liftover/BuildConversion.cpp \
	liftover/Chain.h \
	liftover/Chain.cpp \
	liftover/Converter.h \
	liftover/Converter.cpp \
	liftover/ConverterSQLite.h \
	liftover/ConverterSQLite.cpp \
	liftover/RegionConverter.h \
	liftover/RegionConverter.cpp
	
	
AM_CPPFLAGS=-I$(top_srcdir)/src $(SQLITE_CFLAGS)
#AM_LDFLAGS=-lsqlite3

libknowledge_la_LIBADD=$(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) $(SQLITE_LIBS)