SOCI_HOME=/usr/local
# -- Some libraries and applications require freetype. They will definie this variable
                                                                                                            
ifeq ($(FREETYPE), 1)
ifeq ($(WIN32), 1)
    FT_CPPFLAGS:=-I/c/unx/3rdparty/include -I/c/unx/3rdparty/include/freetype2
    FT_LINK:=-L/c/unx/3rdparty/lib -lfreetype
else
    FT_CPPFLAGS:=$(shell freetype-config --cflags)
    FT_LINK:=$(shell freetype-config --libs)
endif
endif

# -- Some apps depend on wx-widgets
ifeq ($(WX_WIDGETS), 1)
    WX_CPPFLAGS:=$(shell wx-config --cppflags)
    WX_LINK:=$(shell wx-config --libs)
endif


ifeq ($(SOCI), 1)
	  SOCI_CPPFLAGS:=-I$(SOCI_HOME)/soci -I$(SOCI_HOME)/include/soci -I$(SOCI_HOME)/soci/sqlite3 -I$(SOCI_HOME)/include/soci/sqlite3/
	  SOCI_LINKS:=-L$(SOCI_HOME)/lib -lsoci_sqlite3 -lsoci_core -lsqlite3 -ldl -lpthread
endif

EXTERNAL_LINKS+=$(WX_LINK) $(FT_LINK) $(SOCI_LINKS) $(APP_LINKS)
EXT_INCLUDES+=$(WX_CPPFLAGS) $(FT_CPPFLAGS) $(SOCI_CPPFLAGS)
