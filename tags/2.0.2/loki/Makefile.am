SUBDIRS = loadPops

noinst_LTLIBRARIES=liblokiupdate.la

liblokiupdate_la_SOURCES=loki_db.py \
loki-build.py \
loki_source.py \
loki_db.py \
loki_updater.py \
__init__.py \
loaders/__init__.py \
loaders/loki_source_biogrid.py \
loaders/loki_source_dbsnp.py \
loaders/loki_source_entrez.py \
loaders/loki_source_go.py \
loaders/loki_source_kegg.py \
loaders/loki_source_mint.py \
loaders/loki_source_netpath.py \
loaders/loki_source_pfam.py \
loaders/loki_source_pharmgkb.py \
loaders/loki_source_reactome.py \
loaders/loki_source_chainfiles.py \
util/__init__.py \
util/liftOver.py

EXTRA_DIST=loadPops/buildPopulations.py

bin_SCRIPTS=loadPops/buildPopulations.py

if LOKI_INSTALL
VAR_FILES=knowledge.bio
LOKI_FLAGS="-u -v"
else
VAR_FILES=test.bio
LOKI_FLAGS=-t -u -v
endif

LOKI_CMDS=$(srcdir)/loki-build.py -k $(VAR_FILES) $(LOKI_FLAGS); if test "$$?" -ne 0 -a -f "$(VAR_FILES)" ; then mv "$(VAR_FILES)" "$(VAR_FILES).error"; fi

knowledge.bio: liblokiupdate.la .creation
	$(LOKI_CMDS)
	

test.bio: liblokiupdate.la .creation
	$(LOKI_CMDS)
	
.creation:
	date +"%Y.%m.%d" > .creation

DATE_PREFIX:=$(shell cat .creation).

install-data-hook:
	for f in $(VAR_FILES); do \
		mv -f $(DESTDIR)$(datadir)/$$f $(DESTDIR)$(datadir)/$(DATE_PREFIX)$$f; \
		chmod a-wx $(DESTDIR)$(datadir)/$(DATE_PREFIX)$$f; \
		ln -s $(DESTDIR)$(datadir)/$(DATE_PREFIX)$$f $(DESTDIR)$(datadir)/$$f; \
	done

uninstall-hook:
	rm -f $(DESTDIR)$(datadir)/*test.bio $(DESTDIR)$(datadir)/*variations* $(DESTDIR)$(datadir)/*knowledge.bio

nodist_data_DATA=$(VAR_FILES)

CLEANFILES=test.bio knowledge.bio .creation
