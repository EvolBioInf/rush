all: rush

rush: bin/rush
bin/rush:
	make -C rush
	mkdir -p bin
	cp rush/rush bin/
doc: doc/rushDoc.pdf
doc/rushDoc.pdf:
	make -C doc
test: rush
	make test -C rush
clean:
	make clean -C rush
	rm -f bin/rush
	make clean -C doc
