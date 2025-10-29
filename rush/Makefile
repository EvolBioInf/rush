CC=gcc
VERSION = $(shell bash ../scripts/getVersion.sh)
DATE = $(shell bash ../scripts/getDate.sh)

#32 bit version
#when the option VER32 is set, the source will be compiled as the 32 bit version
# preprocessor options: VER32 -> 32-bit version, DEBUG -> debug version, UNIX -> unix version (as opposed to windows version); m64 for unix/mac, exclude for linux!
CFLAGS= -O3 -Wall -Wshadow -pedantic -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -std=c99 -DVER32 -DUNIX -I/opt/local/include/ -L/opt/local/lib -DVERSION="\"$(VERSION)\"" -DDATE="\"$(DATE)\""

CFLAGS= -O3 -Wall -Wshadow -pedantic -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -DVER32 -DVERSION="\"$(VERSION)\"" -DDATE="\"$(DATE)\"" \
-I/opt/local/include/ -I../../Lib/Include -L/opt/local/lib/   -L../../Lib/  #-g #-m64
# The source files, object files, libraries and executable name.
SRCFILES= rush.c interface.c eprintf.c sequenceData.c lcpTree.c shulen.c intervalStack.c stringUtil.c prob.c varSd.c
OBJFILES= rush.o interface.o eprintf.o sequenceData.o stringUtil.o lcpTree.o shulen.o intervalStack.o prob.o varSd.o
LIBS= -lm -lgsl -lgslcblas -ldivsufsort
EXECFILE=rush
DIRECTORY=Rush
# The make rule for the executable
.PHONY : all
all : $(EXECFILE)
$(EXECFILE) : $(OBJFILES)
	$(CC) $(CFLAGS) -o $(EXECFILE) $(OBJFILES) $(LIBS)
$(EXECFILE).o: interface.h eprintf.h
interface.o: interface.h
eprintf.o: eprintf.h

doc:
	cd ../Doc; make clean; make pdf; cd ../$(DIRECTORY)_$(VERSION)

# Other Standard make rules
lint : 
	lint $(SRCFILES) | more
clean:
	rm -f *.o *~
tarfile:
	cd ../Doc; make clean; make pdf; cd ../$(DIRECTORY)_$(VERSION)
	mkdir $(DIRECTORY)_$(VERSION)
	tar cvzfh $(EXECFILE)_$(VERSION).tgz $(DIRECTORY)_$(VERSION)
	mv $(EXECFILE)_$(VERSION).tgz ../
	/bin/rm -r $(DIRECTORY)_$(VERSION)

