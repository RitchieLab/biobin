lib_LTLIBRARIES=libknowledge.la

libknowledge_la_SOURCES=chromosome.cpp \
	def.cpp \
	genegenemodel.cpp \
	genegenemodelarchive.cpp \
	groupmanager.cpp \
	groupmanagerdb.cpp \
	region.cpp \
	regioncontainer.cpp \
	regionmanager.cpp \
	regionmanagerdb.cpp \
	snpdataset.cpp \
	snpsnpmodel.cpp
	
include_HEADERS=chromosome.h \
	genegenemodel.h \
	genegenemodelarchive.h \
	groupmanager.h \
	groupmanagerdb.h \
	region.h \
	regioncontainer.h \
	regionmanager.h \
	regionmanagerdb.h \
	snpdataset.h \
	snpsnpmodel.h 
	
AM_CPPFLAGS=-I..	

libknowledge_la_LIBADD=$(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB)
