base_url=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/

all: ENCODE_human

ENCODE_human: peak_files wgEncodeGencodeBasicV12.gp.gz conversion_table.txt cv.ra
	../../../bin/parse_encode.pl 10000 > ../../../package_data/GeneSets.Homo.sapiens/data/ENCODE.txt
	rm -r peak_files

peak_files: files.txt
	mkdir -p $@
	-cut -f 1 $< | wget -c --quiet --directory-prefix=$@ --base=$(base_url) -i -

wgEncodeGencodeBasicV12.gp.gz:
	wget --quiet http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeGencodeV12/wgEncodeGencodeBasicV12.gp.gz

files.txt:
	wget --quiet $(base_url)files.txt

conversion_table.txt:
	../../../bin/create_conversion_table.R Homo.sapiens ENSEMBLTRANS ENTREZID > $@

cv.ra:
	wget --quiet http://hgdownload.cse.ucsc.edu/goldenPath/encodeDCC/cv.ra

.INTERMEDIATE: files.txt wgEncodeGencodeBasicV12.gp.gz conversion_table.txt cv.ra
