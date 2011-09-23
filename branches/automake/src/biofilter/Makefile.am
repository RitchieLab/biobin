lib_LTLIBRARIES=libbiofilter.la

libbiofilter_la_SOURCES=application.cpp \
	configuration.cpp \
	ldsplineimporter.cpp \
	task.cpp \
	taskgenecoverage.cpp
	
include_HEADERS=application.h \
	configuration.h \
	ldsplineimporter.h \
	task.h \
	taskgenecoverage.h
	
AM_CPPFLAGS=-I..

libbiofilter_la_LIBADD=$(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_REGEX_LIB)
