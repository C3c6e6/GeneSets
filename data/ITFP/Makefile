datasets=human mouse rat
data_files=$(patsubst %, %_tf.txt, $(datasets))
target_files=$(patsubst %, %_tf_regulated_gene.txt, $(datasets))

all: itfp

itfp: $(data_files) $(target_files)
	../../bin/parse_itfp.pl ../../organisms ../../package_data/
	rm *.zip.?

%_tf.txt: %.zip
	unzip -qq $< $@

%_tf_regulated_gene.txt: %.zip
	wget --quiet http://itfp.biosino.org/itfp/system_data/$*/$*.zip
	unzip -qq $*.zip $*_tf_regulated_gene/*
	cat $*_tf_regulated_gene/* > $@
	rm -r $*_tf_regulated_gene

%.zip:
	wget --quiet http://itfp.biosino.org/itfp/system_data/$*/$*.zip

.INTERMEDIATE: $(data_files) $(target_files)

