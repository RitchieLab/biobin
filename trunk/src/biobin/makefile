PROJECT_COMPILER_FLAGS=$(FLAGS_DOPT) $(FLAGS_CV) $(FLAGS_MPI)

PROJECTNAME=biobin

HDRS=binapplication.h \
   binmanager.h \
   configuration.h \
   dataimporter.h \
   genotypestorage.h \
   individual.h \
   main.h \
   taskbincollapse.h \
   taskfilegeneration.h \
   timestamp.l


#ADLDEPS=
#timestamp.l

#DEFINE THE sources
SRCS=binapplication.cpp \
   binmanager.cpp \
   configuration.cpp \
   dataimporter.cpp \
   genotypestorage.cpp \
   main.cpp \
   individual.cpp \
   taskbincollapse.cpp \
   taskfilegeneration.cpp


TEST_SRCS=binmanager.cpp \
   dataimporter.cpp \
   dataset.cpp \
   individual.cpp 

#DEFINE THE ROOTPATH - This is the top level where there will be a make, dist, and others
ROOTPATH=../..

#OVERRIDE the path to sources if necessary
SRCPATHS=..


#list the projects which build into libraries
LOCAL_LIBS=biofilter knowledge utility liftover ldspline vcftools
#ldspline
#appinterface.l

#Uncomment this line for executables
LIBRARY_MAKE=0

#Uncomment this line for test applications
IS_TEST=0
SOCI=1
MAKEPATH=$(ROOTPATH)
EXTERNAL_LINKS=-lboost_filesystem -lboost_system -fprofile-arcs -ftest-coverage -lz -lboost_regex

include $(MAKEPATH)/extlibs.make
include $(MAKEPATH)/make.base


OBJECTONLY: $(OBJS)
	$(SUPRESS)mkdir -p $(HEADERDEST)
	$(SUPRESS)cp $(SRCPATH)/*.h $(HEADERDEST)
	



