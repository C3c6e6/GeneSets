all: PhosphoSitePlus

PhosphoSitePlus: Kinase_Substrate_Dataset
	../../bin/parse_phospho.pl ../../organisms ../../package_data

Kinase_Substrate_Dataset:
	wget --quiet http://www.phosphosite.org/downloads/Kinase_Substrate_Dataset.gz
	gunzip $@.gz

.INTERMEDIATE: Kinase_Substrate_Dataset
