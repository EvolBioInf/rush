all: cli.sh test.sh

cli.sh: intro.org
	bash ../scripts/org2nw intro.org | notangle -Rcli > cli.sh

test.sh: intor.org
	bash ../scriptsorg2nw intro.org | notangle -Rtest.sh > test.sh

test: cli.sh
	bash cli.sh

clean:
	rm -f cli.sh
