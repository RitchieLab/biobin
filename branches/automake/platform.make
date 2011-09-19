
ISAPPLE=$(strip $(shell $(CXX) --version | grep "apple" | wc -l))
ifeq ($(ISAPPLE), 1)
	PLATFORM=APPLE
	ARCH ?= $(shell uname -a | cut -d " " -f 15)
	BITCOUNT ?= $(shell getconf LONG_BIT)
	LARGE_FILESYSTEM=-D_FILE_OFFSET_BITS=64
	PLAT_INCLUDES=-I/usr/X11/include
	PLAT_LINKS=-L/usr/X11/lib
	EXE_PLAT=-OSX
endif


ifeq ($(shell uname), Linux)
	ARCH ?= $(shell uname -a | cut -d " " -f 12)
	BITCOUNT ?= $(shell getconf LONG_BIT)
	ifneq ($(BITCOUNT), 32)
		EXE_PLAT=$(BITCOUNT)
		LARGE_FILESYSTEM=-D_FILE_OFFSET_BITS=64
	endif
	ifeq ($(shell $(CXX) --version | grep "Red Hat" | wc -l), 1)
		PLATFORM=RH
	else
		PLATFORM=UX
	endif
endif
ifdef ($(ARCH))
	MARCH=-march=$(ARCH)
endif

ifdef WIN32
	#GNUCC=/home/torstenson/MinGW/bin/i686-pc-mingw32-cpp
	GNUCC=g++
	ARCH=i686
	BITCOUNT=32
	EXE_PLAT=
	PLATFORM=win
	PLAT_INCLUDES=-I/c/unx/3rdparty/include -I/usr/local/include
else
	GNUCC=g++
endif

ifdef ($BITCOUNT)
	WORDSIZE=-m$(BITCOUNT)
endif

