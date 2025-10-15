date = $(shell git log | grep -m 1 Date | sed -r 's/Date: +[A-Z][a-z]+ ([A-Z][a-z]+) ([0-9]+) [^ ]+ ([0-9]+) .+/\2_\1_\3/')
version = $(shell git describe)

all: rushDoc.pdf

rushDoc.pdf: rushDoc.tex header.tex intro.tex
	echo $(date) | tr '_' ' ' > date.txt
	echo $(version) | tr '-' ' ' | awk '{printf "%s", $$1; if ($$2) printf "-%s", $$2; printf "\n"}' > version.txt
	latex rushDoc
	bibtex rushDoc
	latex rushDoc
	latex rushDoc
	dvips -o -q rushDoc
	ps2pdf -dALLOWPSTRANSPARENCY rushDoc.ps
intro.tex: ../intro/intro.org
	bash ../scripts/org2nw ../intro/intro.org | noweave -n -x > intro.tex
publish: rushDoc.pdf
	if [ -d ~/WinHome/ownCloud ]; then \
		cp rushDoc.pdf ~/WinHome/ownCloud/docs/; \
	fi
clean:
	rm -f rushDoc.aux rushDoc.bbl rushDoc.bib rushDoc.blg rushDoc.d rushDoc.dvi \
	rushDoc.log rushDoc.pdf rushDoc.ps
