noinst_LTLIBRARIES=libbiofilter.la

libbiofilter_la_SOURCES=application.cpp \
	configuration.cpp \
	ldsplineimporter.cpp \
	task.cpp \
	taskgenecoverage.cpp \
	application.h \
	configuration.h \
	ldsplineimporter.h \
	task.h \
	taskgenecoverage.h \
	taskcoveragereport.h \
taskgenebased.h \
taskgenegenemodelreport.cpp \
taskgenegenemodelreport.h \
taskgenereport.h \
taskmarkerinfo.h \
tasksnpgenemap.h \
tasksnpreport.h \
tasksnpsnpmodelarchive.h
	
AM_CPPFLAGS=-I$(top_srcdir)/src -DDATA_DIR='"$(datadir)"'

libbiofilter_la_LIBADD=$(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_REGEX_LIB)
