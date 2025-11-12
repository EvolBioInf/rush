CC=gcc
VERSION = $(shell bash ../scripts/getVersion.sh)
DATE = $(shell bash ../scripts/getDate.sh)

CFLAGS= -O3 -Wall -Wshadow -pedantic -DVERSION="\"$(VERSION)\"" -DDATE="\"$(DATE)\""

# The source files, object files, libraries and executable name.
SRCFILES= rush.c interface.c eprintf.c sequenceData.c lcpTree.c shulen.c intervalStack.c stringUtil.c prob.c varSd.c
OBJFILES= rush.o interface.o eprintf.o sequenceData.o stringUtil.o lcpTree.o shulen.o intervalStack.o prob.o varSd.o
LIBS= -lm -lgsl -lgslcblas -ldivsufsort

EXECFILE=rush
# The make rule for the executable
.PHONY : all
all : $(EXECFILE)
$(EXECFILE) : $(OBJFILES)
	$(CC) $(CFLAGS) -o $(EXECFILE) $(OBJFILES) $(LIBS)
$(EXECFILE).o: interface.h eprintf.h
interface.o: interface.h
eprintf.o: eprintf.h

test: $(EXECFILE)
	bash makeTest.sh
clean:
	rm -f *.o *~
