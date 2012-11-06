ACLOCAL_AMFLAGS= -I ./make/m4

SUBDIRS=src/vcftools src/knowledge src/biobin updater

EXTRA_DIST=CREATED examples

CLEANFILES=CREATED

dist-hook:
	chmod u+w $(distdir)/CREATED
	echo "Distribution date: `date`" > CREATED

# This is here to prevent installation of LOKI during a distcheck
DISTCHECK_CONFIGURE_FLAGS=--disable-loki
