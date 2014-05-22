organisms=$(shell awk 'gsub("[ -]", ".")' organisms | awk '{gsub(/(sub)?str\./, ""); print}' | awk 'gsub(/\.+/, ".")')
targets=$(patsubst %, %.built, $(organisms))

all: $(targets)

$(targets): %.built: package_data/GeneSets.%/checksums
	bin/create_package.pl $(<D)
	mv package_data/GeneSets.$*_*.tar.gz ./
	touch $@

package_data/GeneSets.%/checksums: datafiles
	bin/checksums.py $(@D)
	
package_data: organisms
	bin/create_package_dirs

datafiles: package_data
	$(MAKE) -C data/
	touch datafiles

.INTERMEDIATE: datafiles

clean:
	rm -fr package_data *.built *.tar.gz
