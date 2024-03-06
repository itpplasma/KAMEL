# load variables defined in top level Makefile.in
include ../../Makefile.in

#targets:
all: check test

check:
ifeq ($(suf),)
	$(error Error: a language for the test sources is not correctly specified: $(lang). Call: make lang=c++ or make lang=c or make lang=fortran)
endif

test:   test.o
	$(COMP) -o test.$(suf).x $< $(LDFLAGS)
	rm -f test.o

clean:
	rm -f *.o *.x *.dump

test.o:
	$(COMP) $(COMPOPTS) -c test.$(suf)
