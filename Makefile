all: rush

rush: bin/rush
bin/rush:
	make -C rush
	cp rush/rush bin
doc: doc/rushDoc.pdf
doc/rushDoc.pdf:
	make -C doc
clean:
	make clean -C rush
	rm bin/rush
	make clean -C doc
