all: bin/rush doc

bin/rush:
	make -C rush
	cp rush/rush bin
doc:
	make -C doc
clean:
	make clean -C rush
	make clean -C doc
