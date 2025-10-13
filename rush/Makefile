CC=gcc
CFLAGS= -O3 -Wall -Wshadow -pedantic -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -DVER32 \
-I/opt/local/include/ -I../../Lib/Include -L/opt/local/lib/   -L../../Lib/  #-g #-m64
# The source files, object files, libraries and executable name.
SRCFILES= rush.c interface.c eprintf.c sequenceData.c lcpTree.c shulen.c intervalStack.c stringUtil.c minimize.c prob.c varSd.c
OBJFILES= rush.o interface.o eprintf.o sequenceData.o stringUtil.o lcpTree.o shulen.o intervalStack.o \
	DeepShallow64/ds_ssort.a DeepShallow64/bwtlcp.a minimize.o prob.o varSd.o
LIBS= -lm DeepShallow64/ds_ssort.a DeepShallow64/bwtlcp.a -lgsl -lgslcblas
EXECFILE=rush
DIRECTORY=Rush
VERSION=1.4
# The make rule for the executable
.PHONY : all
all : $(EXECFILE) DeepShallow64/ds_ssort.a DeepShallow64/bwtlcp.a
$(EXECFILE) : $(OBJFILES)
	$(CC) $(CFLAGS) -o $(EXECFILE) $(OBJFILES) $(LIBS)
$(EXECFILE).o: interface.h eprintf.h
interface.o: interface.h
eprintf.o: eprintf.h
DeepShallow64/ds_ssort.a:
	cd DeepShallow64; make
DeepShallow64/bwtlcp.a:
	cd DeepShallow64; make

doc:
	cd ../Doc; make clean; make pdf; cd ../$(DIRECTORY)_$(VERSION)

# Other Standard make rules
lint : 
	lint $(SRCFILES) | more
clean:
	rm -f *.o *~; cd DeepShallow64; make clean
tarfile:
	cd ../Doc; make clean; make pdf; cd ../$(DIRECTORY)_$(VERSION)
	mkdir $(DIRECTORY)_$(VERSION)
	cd DeepShallow64; make clean;
	cp -rf DeepShallow64 Scripts $(SRCFILES) *.h  README COPYING Makefile ../Doc/rushDoc.pdf $(DIRECTORY)_$(VERSION)
	tar cvzfh $(EXECFILE)_$(VERSION).tgz $(DIRECTORY)_$(VERSION)
	mv $(EXECFILE)_$(VERSION).tgz ../
	/bin/rm -r $(DIRECTORY)_$(VERSION)

