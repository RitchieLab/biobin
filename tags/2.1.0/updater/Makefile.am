SUBDIRS = loki/loadPops

noinst_LTLIBRARIES=liblokiupdate.la

liblokiupdate_la_SOURCES=\
loki-build.py \
loki/loki_source.py \
loki/loki_db.py \
loki/loki_updater.py \
loki/__init__.py \
loki/loaders/__init__.py \
loki/loaders/loki_source_biogrid.py \
loki/loaders/loki_source_dbsnp.py \
loki/loaders/loki_source_entrez.py \
loki/loaders/loki_source_go.py \
loki/loaders/loki_source_kegg.py \
loki/loaders/loki_source_mint.py \
loki/loaders/loki_source_netpath.py \
loki/loaders/loki_source_pfam.py \
loki/loaders/loki_source_pharmgkb.py \
loki/loaders/loki_source_reactome.py \
loki/loaders/loki_source_chainfiles.py \
loki/loaders/loki_source_gwas.py \
loki/loaders/loki_source_ucsc_ecr.py \
loki/loaders/loki_source_oreganno.py \
loki/loaders/test/__init__.py \
loki/loaders/test/loki_source_genes.py \
loki/loaders/test/loki_source_light.py \
loki/loaders/test/loki_source_paint.py \
loki/loaders/test/loki_source_snps.py \
loki/loaders/test/loki_source_spectrum.py \
loki/util/__init__.py \
loki/util/liftOver.py

EXTRA_DIST=loki/loadPops/buildPopulations.py

bin_SCRIPTS=loki/loadPops/buildPopulations.py

if LOKI_INSTALL
VAR_FILES=knowledge.bio
LOKI_FLAGS=-u -v
else
VAR_FILES=test.bio
LOKI_FLAGS=-t -u -v
endif

LOKI_CMDS=$(srcdir)/loki-build.py -k $(VAR_FILES) $(LOKI_FLAGS); if test "${?}" -ne 0 -a -f "$(VAR_FILES)" ; then mv "$(VAR_FILES)" "$(VAR_FILES).error"; fi

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
