# List of subdirectories
SUBDIRS = KiLCA KIM QL-Balance

# Default target: build all
all: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

clean:
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir clean; done

.PHONY: all clean $(SUBDIRS)
