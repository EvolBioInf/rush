all: test.sh

test.sh: intro.org
	bash ../scripts/org2nw intro.org | notangle -Rtest.sh > test.sh

test: test.sh
	bash test.sh

clean:
	rm -f test.sh
