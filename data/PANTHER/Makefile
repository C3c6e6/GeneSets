all: PANTHER

PANTHER: BioPAX.tar.gz ../BioSystems/names.dmp
	tar -xzvf $<
	Rscript ../../bin/parse_BioPax3.R BioPAX $@ ../../organisms ../BioSystems/names.dmp ../../package_data
	rm -r BioPAX

../BioSystems/names.dmp:
	$(MAKE) -C $(@D) $(@F)

BioPAX.tar.gz: 
	wget --quiet ftp://ftp.pantherdb.org//pathway/current_release/BioPAX.tar.gz

.INTERMEDIATE: BioPAX.tar.gz ../BioSystems/names.dmp
