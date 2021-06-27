# Makefile for the mp documentation


# User Guide and Reference manual
mp: mp.tex
	pdflatex mp
	makeindex mp.idx
	bibtex mp

#	mv mp.pdf $(WEBPAGE)/documentation


# mp reference guide
mpRef: mpRef.tex
	pdflatex mpRef
	makeindex mpRef.idx
	bibtex mpRef

# Notes on solving multi-domain problems
multiDomain: multiDomain.tex
	latex multiDomain
	makeindex multiDomain.idx
	bibtex multiDomain
	latex multiDomain
	dvips -P pdf -G0 -o multiDomain.ps multiDomain.dvi
	ps2pdf multiDomain.ps


# Notes on solving fluid solid-mechanics problems
fsm: fluidStructure/fsm.tex
	cd fluidStructure; make fsm


fig:
	latex fig.tex
	dvips -P pdf -G0 -o fig.ps fig.dvi
	ps2pdf fig.ps

clean: 
	-rm -f *.aux *.bbl *.blg *.ind *.ilg *.log *.idx *.toc

