# config, can be set via env variable ('CONFIG') before running make
CONFIG ?= Release

# List of subdirectories
SUBDIRS := KiLCA KIM QL-Balance

# Default target: build all
all: $(SUBDIRS)

$(SUBDIRS):
	export CONFIG=$(CONFIG)
	$(MAKE) -C $@

clean:
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir clean; done

.PHONY: all clean $(SUBDIRS)
