PROJECT_COMPILER_FLAGS=$(FLAGS_DOPT) $(FLAGS_CV) $(FLAGS_MPI)

PROJECTNAME=utility

HDRS=filetools.h \
	lineparser.h \
	locus.h \
	strings.h
#DEFINE THE sources
SRCS=filetools.cpp \
	lineparser.cpp \
	locus.cpp \
	strings.cpp 

#DEFINE THE ROOTPATH - This is the top level where there will be a make, dist, and others
ROOTPATH=../..

#OVERRIDE the path to sources if necessary
SRCPATHS=..

#indicate what the executable will be named
#EXENAME=ese

#list the projects which build into libraries
#LOCAL_LIBS=genetics utility

#Uncomment this line for executables
LIBRARY_MAKE=1

#Uncomment this line for test applications
IS_TEST=0



MAKEPATH=$(ROOTPATH)
include $(MAKEPATH)/make.base


OBJECTONLY: $(OBJS)
	$(SUPRESS)mkdir -p $(HEADERDEST)
	$(SUPRESS)cp $(SRCPATH)/*.h $(HEADERDEST)
	



