noinst_LTLIBRARIES=libliftover.la

libliftover_la_SOURCES=chain.cpp \
	converter.cpp \
	converterdb.cpp \
	region.cpp \
	chain.h \
	converter.h \
	converterdb.h \
	region.h \
	snp.h
		
AM_CPPFLAGS=-I$(top_srcdir)/src

libliftover_la_LIBADD=$(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB)
