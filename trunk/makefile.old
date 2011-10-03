BINPATH=./bin
ROOTPATH=.

include $(ROOTPATH)/platform.make
include $(ROOTPATH)/color.make

ifeq ($(INC_VERSION), 1)
	INCREMENT=INC_BUILD=1
endif 

ifeq ($(DEBUG), 1)
	LIBSUFFIX=d
	RELEASEMODE=RELEASE=0
	RELEASE_DESC="Debug Mode"
	INCREMENT=
else
	LIBSUFFIX=
	RELEASEMODE=RELEASE=1
	RELEASE_DESC="Release Mode"
endif

ifeq ($(USE_DOPT), 1)
	OPTION_DOPT=USE_DOPT=1
endif

ifeq ($(ENABLE_CV), 0)
	OPTION_CV=
	CV_SUFF=
else
	OPTION_CV=ENABLE_CV=1
	CV_SUFF="_cv"
endif
GNUCC=g++

ifdef COMPAT
	CPP=$(GNUCC)$(COMPAT)
else
	CPP=$(GNUCC)
endif

ifeq ($(LOCAL_RELEASE), 1)
	CPPFLAGS~=LOCAL_RELEASE=1
endif


ifeq ($(USE_MPI), 1)
	FLAGS_MPI=-D USE_MPI
	MPI_LINK=$(MPI_LIBS)
	EXEC_PREFIX=p
	MPIPATH=parallel
	MPI_DESC="w/ MPI"
else
	MPI_LINK=
	EXEC_PREFIX=
	MPIPATH=serial
	MPI_DESC="w/o MPI"
endif

include $(ROOTPATH)/platform.make
#LIBPATH=$(ROOTPATH)/bin/lib

LIBPATH=$(ROOTPATH)/bin/lib/$(MPIPATH)/$(PLATFORM)-$(BITCOUNT)bit/$(CPP)

GCC_VERSION=$(shell $(CXX) --version | grep gcc | cut -f3 -d\ | cut -f1-2 -d. )


ifdef VAMPIRE
    MPI_LIBS= -lmpich -lpthread
    MPI_INC=/usr/local/mpich-ether/p4/absoft_gcc/include/
    CL=mpiCC
	VMP_SUFF="_vmp"
else
	MPI_LIBS=-lmpichcxx -lmpich -lpthread
    MPI_INC=
    CL=$(CPP)
	VMP_SUFF=
endif

ifeq ($(VERBOSE), 1)
   PROJECT_OPTIONS=$(OPTION_DOPT) $(OPTION_CV) $(RELEASE_MODE) $(INCREMENT) VERBOSE=1
else
   PROJECT_OPTIONS=-s $(OPTION_DOPT) $(OPTION_CV) $(RELEASE_MODE) $(INCREMENT)
endif

LOCAL_LIBS=knowledge utility ldspline liftover biofilter vcftools
#ldspline pngwriter utility genetics simulation simla simpen ga random


DEPENDANTLIBS=$(patsubst %, $(LIBPATH)/lib%$(EXEC_PREFIX)$(LIBSUFFIX).a, $(LOCAL_LIBS))

all: biofilterapp

.PHONY: clean test 


utility: $(LIBPATH)/libutility$(LIBSUFFIX).a
random: $(LIBPATH)/librandom$(LIBSUFFIX).a
genomesim: $(LIBPATH)/libgenomesim$(LIBSUFFIX).a
liftover: $(LIBPATH)/libliftover$(LIBSUFFIX).a
pngwriter: $(LIBPATH)/libpngwriter$(LIBSUFFIX).a
vcftools: $(LIBPATH)/vcftools$(LIBSUFFIX).a
ldspline: $(BINPATH)/ldspline$(LIBSUFFIX)
biobin: $(BINPATH)/biobin$(LIBSUFFIX)
biofilter: $(LIBPATH)/libbiofilter$(EXEC_PREFIX)$(LIBSUFFIX).a
biofilterapp: $(BINPATH)/biofilter$(LIBSUFFIX)
paris: $(BINPATH)/paris$(LIBSUFFIX)
knowledge: $(LIBPATH)/libknowledge$(EXEC_PREFIX)$(LIBSUFFIX).a
paris-test: $(BINPATH)/paris-test$(LIBSUFFIX)
biofilter-test: $(BINPATH)/test-biofilter$(LIBSUFFIX)
utility-test:  $(BINPATH)/test-utility$(LIBSUFFIX)
knowledge-test: $(BINPATH)/test-knowledge$(LIBSUFFIX)
liftover-test: $(BINPATH)/test-liftover$(LIBSUFFIX)


release:
	make COMPAT=32 

debug: 
	make DEBUG=1 

cleantests:
	@echo "Removing test applications"
	@rm -f $(BINPATH)/test-utility$(EXE_PLAT)$(LIBSUFFIX)
	@rm -f $(BINPATH)/test-knowledge$(EXE_PLAT)$(LIBSUFFIX)
	@rm -f $(BINPATH)/test-biofilter$(EXE_PLAT)$(LIBSUFFIX)
clean:
	@$(MAKE) -C $(ROOTPATH)/src/utility $(RELEASEMODE) clean
	@$(MAKE) -C $(ROOTPATH)/src/liftover $(RELEASEMODE) clean
	@$(MAKE) -C $(ROOTPATH)/src/random $(RELEASEMODE) clean
	@$(MAKE) -C $(ROOTPATH)/src/vcftools $(RELEASEMODE) clean
	@$(MAKE) -C $(ROOTPATH)/src/biofilter $(RELEASEMODE) clean
	@$(MAKE) -C $(ROOTPATH)/src/knowledge $(RELEASEMODE) clean
	@$(MAKE) -C $(ROOTPATH)/src/vcftools $(RELEASEMODE) clean

cleantmps:
	@$(MAKE) -C $(ROOTPATH)/src/utility $(RELEASEMODE) cleantmps
	@$(MAKE) -C $(ROOTPATH)/src/liftover $(RELEASEMODE) cleantmps
	@$(MAKE) -C $(ROOTPATH)/src/random $(RELEASEMODE) cleantmps
	@$(MAKE) -C $(ROOTPATH)/src/vcftools $(RELEASEMODE) cleantmps
	@$(MAKE) -C $(ROOTPATH)/src/biofilter $(RELEASEMODE) cleantmps
	@$(MAKE) -C $(ROOTPATH)/src/knowledge $(RELEASEMODE) cleantmps
	@$(MAKE) -C $(ROOTPATH)/src/vcftools $(RELEASEMODE) cleantmps

tests: cleantests utility-test knowledge-test biofilter-test liftover-test

$(BINPATH)/biobin$(LIBSUFFIX): $(DEPENDANTLIBS)
	@echo $(MAKE) -C $(ROOTPATH)/src/biobin $(RELEASEMODE) $(PROJECT_OPTIONS)
	@$(MAKE) -C $(ROOTPATH)/src/biobin $(RELEASEMODE) $(PROJECT_OPTIONS)

$(BINPATH)/biofilter$(LIBSUFFIX): $(DEPENDANTLIBS)
	@echo "LIBS: $(DEPENDANTLIBS)"
	@echo $(MAKE) -C $(ROOTPATH)/src/biobin $(RELEASEMODE) $(PROJECT_OPTIONS)
	@$(MAKE) -C $(ROOTPATH)/src/biofilter -f make.biofilter $(RELEASEMODE) $(PROJECT_OPTIONS)
$(BINPATH)/test-biofilter$(LIBSUFFIX): $(DEPENDANTLIBS)
	@$(MAKE) -C $(ROOTPATH)/src/biofilter -f make.test $(RELEASEMODE) $(PROJECT_OPTIONS) IS_TEST=1
$(BINPATH)/ldspline$(LIBSUFFIX): $(DEPENDENTLIBS) $(LIBPATH)/libldspline$(EXEC_PREFIX)$(LIBSUFFIX).a
	@$(MAKE) -C $(ROOTPATH)/src/ldspline -f main.make $(RELEASEMODE) $(PROJECT_OPTIONS)
$(BINPATH)/test-utility$(LIBSUFFIX):
	@$(MAKE) -C $(ROOTPATH)/src/utility -f make.test $(RELEASEMODE) $(PROJECT_OPTIONS) IS_TEST=1
$(BINPATH)/test-knowledge$(LIBSUFFIX):
	@$(MAKE) -C $(ROOTPATH)/src/knowledge -f make.test $(RELEASEMODE) $(PROJECT_OPTIONS) IS_TEST=1
$(BINPATH)/test-liftover$(LIBSUFFIX):
	@$(MAKE) -C $(ROOTPATH)/src/liftover -f make.test $(RELEASEMODE) $(PROJECT_OPTIONS) IS_TEST=1
$(LIBPATH)/libvcftools$(EXEC_PREFIX)$(LIBSUFFIX).a:
	$(TERM_CYAN)
	@echo "Building $@ $^"
	$(TERM_RESET)
	@$(MAKE) -C $(ROOTPATH)/src/vcftools $(RELEASEMODE) $(PROJECT_OPTIONS)
	
$(LIBPATH)/libbiofilter$(EXEC_PREFIX)$(LIBSUFFIX).a:
	$(TERM_CYAN)
	@echo "Building $@ $^"
	$(TERM_RESET)
	@$(MAKE) -C $(ROOTPATH)/src/biofilter $(RELEASEMODE) $(PROJECT_OPTIONS)

$(LIBPATH)/libldspline$(EXEC_PREFIX)$(LIBSUFFIX).a:
	@$(MAKE) -C $(ROOTPATH)/src/ldspline $(RELEASEMODE) $(PROJECT_OPTIONS)

$(LIBPATH)/libgenetics$(EXEC_PREFIX)$(LIBSUFFIX).a:
	$(TERM_CYAN)
	@echo "Building $@ $^"
	$(TERM_RESET)
	@$(MAKE) -C $(ROOTPATH)/src/genetics $(RELEASEMODE) $(PROJECT_OPTIONS)
	
$(LIBPATH)/libknowledge$(EXEC_PREFIX)$(LIBSUFFIX).a:
	$(TERM_CYAN)
	@echo "Building $@ $^"
	$(TERM_RESET)
	@$(MAKE) -C $(ROOTPATH)/src/knowledge $(RELEASEMODE) $(PROJECT_OPTIONS)

$(LIBPATH)/libliftover$(LIBSUFFIX).a:
	$(TERM_CYAN)
	@echo "Building $@ $^"
	$(TERM_RESET)
	@$(MAKE) -C $(ROOTPATH)/src/liftover $(RELEASEMODE) $(PROJECT_OPTIONS)

$(LIBPATH)/librandom$(EXEC_PREFIX)$(LIBSUFFIX).a:
	$(TERM_CYAN)
	@echo "Building $@ $^"
	$(TERM_RESET)
	@$(MAKE) -C $(ROOTPATH)/src/random $(RELEASEMODE) $(PROJECT_OPTIONS)
	
$(LIBPATH)/libpngwriter$(EXEC_PREFIX)$(LIBSUFFIX).a:
	$(TERM_CYAN)
	@echo "Building $@ $^"
	$(TERM_RESET)
	@$(MAKE) -C $(ROOTPATH)/src/pngwriter $(RELEASEMODE) $(PROJECT_OPTIONS)

$(LIBPATH)/libsimpen$(EXEC_PREFIX)$(LIBSUFFIX).a:
	$(TERM_CYAN)
	@echo "Building $@ $^"
	$(TERM_RESET)
	@$(MAKE) -C $(ROOTPATH)/src/simpen $(RELEASEMODE) $(PROJECT_OPTIONS)

$(LIBPATH)/libsimla$(EXEC_PREFIX)$(LIBSUFFIX).a:
	$(TERM_CYAN)
	@echo "Building $@ $^"
	$(TERM_RESET)
	@$(MAKE) -C $(ROOTPATH)/src/simla $(RELEASEMODE) $(PROJECT_OPTIONS)

$(LIBPATH)/libga$(EXEC_PREFIX)$(LIBSUFFIX).a:
	$(TERM_CYAN)
	@echo "Building $@ $^"
	$(TERM_RESET)
	@$(MAKE) -C $(ROOTPATH)/src/ga $(RELEASEMODE) $(PROJECT_OPTIONS)

$(LIBPATH)/libsimulation$(EXEC_PREFIX)$(LIBSUFFIX).a:
	$(TERM_CYAN)
	@echo "Building $@ $^"
	$(TERM_RESET)
	@$(MAKE) -C $(ROOTPATH)/src/simulation $(RELEASEMODE) $(PROJECT_OPTIONS)

$(LIBPATH)/libutility$(EXEC_PREFIX)$(LIBSUFFIX).a:
	$(TERM_CYAN)
	@echo "Building $@ $^"
	$(TERM_RESET)
	@$(MAKE) -C $(ROOTPATH)/src/utility $(RELEASEMODE) $(PROJECT_OPTIONS)


$(BINPATH)/test/mdrpdt.test: $(patsubst %, $(LIBPATH)/lib%$(LIBSUFFIX).a, $(LOCAL_LIBS))
	@$(MAKE) -C $(ROOTPATH)/test $(RELEASEMODE) $(PROJECT_OPTIONS)
