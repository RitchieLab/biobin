PROJECT_COMPILER_FLAGS=$(FLAGS_DOPT) $(FLAGS_CV) $(FLAGS_MPI)

PROJECTNAME=ldspline



#DEFINE THE sources
SRCS=ldspline.cpp \
	locuslookup.cpp \
	snpspline.cpp

#DEFINE THE ROOTPATH - This is the top level where there will be a make, dist, and others
ROOTPATH=../..

#OVERRIDE the path to sources if necessary
SRCPATHS=..


#list the projects which build into libraries
LOCAL_LIBS=utility
#appinterface.l

#Uncomment this line for executables
#LIBRARY_MAKE=0

#Uncomment this line for test applications
IS_TEST=0

#appinterface.l: 
	#inc_version.sh timestamp.h
LIBRARY_MAKE=1

MAKEPATH=$(ROOTPATH)
FREETYPE=1
#EXTERNAL_LINKS=-L/opt/local/lib/mysql5/mysql -lmysqlclient
#EXT_INCLUDES=-I/opt/local/include/mysql5/mysql
include $(MAKEPATH)/extlibs.make
include $(MAKEPATH)/make.base



OBJECTONLY: $(OBJS)
	$(SUPRESS)mkdir -p $(HEADERDEST)
	$(SUPRESS)cp $(SRCPATH)/*.h $(HEADERDEST)
	



