organisms=$(shell awk 'gsub(" ", ".")' ../../organisms)

all: $(organisms)

$(organisms): %:
	../../bin/GO.R $* ../../package_data
