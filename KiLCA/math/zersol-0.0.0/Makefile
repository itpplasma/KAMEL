# load variables defined in top level Makefile.in
include Makefile.in

#targets:
all: check bindings

check:
ifeq ($(suf),)
	$(error Error: a language for the bindings is not correctly specified: $(lang). Call: make lang=c, make lang=fortran)
endif

bindings: src/zersol$(suf).o
	$(AR) lib/libzersol$(suf).a src/zersol$(suf).o
	rm -f src/*.o

clean:
	rm -f src/*.o lib/*.a
