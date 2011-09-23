lib_LTLIBRARIES=libliftover.la

libliftover_la_SOURCES=chain.cpp \
	converter.cpp \
	converterdb.cpp \
	region.cpp
	
include_HEADERS=chain.h \
	converter.h \
	converterdb.h \
	region.h 
		
AM_CPPFLAGS=-I..

libliftover_la_LIBADD=$(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB)
