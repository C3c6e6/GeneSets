all: biosystems

base_url=ftp://ftp.ncbi.nih.gov/pub/biosystems/CURRENT/

biosystems: biosystems_gene.gz bsid2info.gz names.dmp gene_info.gz ../../bin/parse_biosystems.pl ../../organisms
	../../bin/parse_biosystems.pl ../../organisms
        
biosystems_gene.gz bsid2info.gz:
	wget --quiet $(base_url)$@

names.dmp: taxdump.tar.gz
	tar -xzf $< $@

taxdump.tar.gz:
	wget --quiet ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

gene_info.gz:
	wget --quiet ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz

.INTERMEDIATE: biosystems_gene.gz bsid2info.gz names.dmp taxdump.tar.gz gene_info.gz
