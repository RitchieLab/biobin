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
	Information.cpp \
	InformationSQLite.h \
	InformationSQLite.cpp \
	Configuration.h \
	Configuration.cpp \
	liftover/Chain.h \
	liftover/Chain.cpp \
	liftover/Converter.h \
	liftover/Converter.cpp \
	liftover/ConverterSQLite.h \
	liftover/ConverterSQLite.cpp \
	liftover/Segment.h \
	any_iterator.hpp \
	detail/any_iterator_abstract_base.hpp \
	detail/any_iterator_metafunctions.hpp \
	detail/any_iterator_wrapper.hpp
	
	
AM_CPPFLAGS=-I$(top_srcdir)/src $(SQLITE_CFLAGS) $(BOOST_CPPFLAGS)
AM_LDFLAGS=$(BOOST_LDFLAGS) 

libknowledge_la_LIBADD=$(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) $(SQLITE_LIBS)
