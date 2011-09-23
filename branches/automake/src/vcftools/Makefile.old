PROJECT_COMPILER_FLAGS=$(FLAGS_DOPT) $(FLAGS_CV) $(FLAGS_MPI)

PROJECTNAME=vcftools

HDRS=dgeev.h vcf_file.h vcf_entry.h output_log.h


#ADLDEPS=
#timestamp.l

#DEFINE THE sources
SRCS=vcf_file.cpp vcf_entry.cpp \
	vcf_entry_getters.cpp vcf_entry_setters.cpp \
	vcf_file_filters.cpp vcf_file_output.cpp \
	vcf_file_format_convert.cpp \
        vcf_file_diff.cpp output_log.cpp

TESTs_SRCS=

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
APP_LINKS=-lgz
include $(MAKEPATH)/extlibs.make
include $(MAKEPATH)/make.base


OBJECTONLY: $(OBJS)
	$(SUPRESS)mkdir -p $(HEADERDEST)
	$(SUPRESS)cp $(SRCPATH)/*.h $(HEADERDEST)
	



