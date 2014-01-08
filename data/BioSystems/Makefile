all: biosystems

base_url=ftp://ftp.ncbi.nih.gov/pub/biosystems/CURRENT/

biosystems: biosystems_gene_all.gz bsid2info.gz names.dmp gene_info.gz go_term_domains.txt ../../bin/parse_biosystems.pl ../../organisms
	../../bin/parse_biosystems.pl ../../organisms
        
biosystems_gene_all.gz bsid2info.gz:
	wget --quiet $(base_url)$@

names.dmp: taxdump.tar.gz
	tar -xzf $< $@

taxdump.tar.gz:
	wget --quiet ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

gene_info.gz:
	wget --quiet ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz

go_term_domains.txt: ../../bin/go_term_domains.R
	../../bin/go_term_domains.R

.INTERMEDIATE: biosystems_gene_all.gz bsid2info.gz names.dmp taxdump.tar.gz gene_info.gz go_term_domains.txt
