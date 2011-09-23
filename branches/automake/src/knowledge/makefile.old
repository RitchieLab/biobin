PROJECT_COMPILER_FLAGS=$(FLAGS_DOPT) $(FLAGS_CV) $(FLAGS_MPI)

PROJECTNAME=knowledge
DEADHRD=snptranslator.h
HDRS=chromosome.h \
	genegenemodel.h \
	genegenemodelarchive.h \
	groupmanager.h \
	groupmanagerdb.h \
	region.h \
	regioncontainer.h \
	regionmanager.h \
	regionmanagerdb.h \
	snpdataset.h \
	snpsnpmodel.h \
	timestamp.l


#ADLDEPS=
#timestamp.l

#DEFINE THE sources
SRCS=chromosome.cpp \
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

TESTs_SRCS=chromosome.cpp \
	def.cpp \
	genegenemodel.cpp \
	genegenemodelarchive.cpp \
	groupmanager.cpp \
	groupmanagerdb.cpp \
	region.cpp \
	regionmanager.cpp \
	regionmanagerdb.cpp \
	testmain.cpp \
	snpdataset.cpp \
	snpsnpmodel.cpp

#DEFINE THE ROOTPATH - This is the top level where there will be a make, dist, and others
ROOTPATH=../..

#OVERRIDE the path to sources if necessary
SRCPATHS=..


#list the projects which build into libraries
#LOCAL_LIBS=utility ldspline
#appinterface.l

#Uncomment this line for executables
LIBRARY_MAKE=1

#Uncomment this line for test applications
IS_TEST=0
SOCI=1
MAKEPATH=$(ROOTPATH)
EXTERNAL_LINKS=
include $(MAKEPATH)/extlibs.make
include $(MAKEPATH)/make.base


OBJECTONLY: $(OBJS)
	$(SUPRESS)mkdir -p $(HEADERDEST)
	$(SUPRESS)cp $(SRCPATH)/*.h $(HEADERDEST)
	



