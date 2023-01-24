# Makefile for the Cgcssi documentation


all = user ref

# User Guide
# obFiles = $(OverBlown)/OverBlown.C $(OverBlown)/setParameters.C $(OverBlown)/getInitialConditions.C \
#           $(OverBlown)/OB_CompositeGridSolver.C $(OverBlown)/OB_Parameters.C
user: CgcssiUserGuide.tex $(CgcssiFiles)
	extract.p $(CgcssiFiles)
	pdflatex CgcssiUserGuide
	bibtex CgcssiUserGuide
	makeindex CgcssiUserGuide.idx

#	mv CgcssiUserGuide.pdf $(WEBPAGE)/documentation
userOld: CgcssiUserGuide.tex $(CgcssiFiles)
	extract.p $(CgcssiFiles)
	latex CgcssiUserGuide
	makeindex CgcssiUserGuide.idx
	dvips -P pdf -G0 -o CgcssiUserGuide.ps CgcssiUserGuide.dvi
	ps2pdf CgcssiUserGuide.ps
#	mv CgcssiUserGuide.pdf $(WEBPAGE)/documentation

# Reference Guide
ref: CgcssiRef.tex
	pdflatex CgcssiRef
	bibtex CgcssiRef
	makeindex CgcssiRef.idx

refOlder: CgcssiRef.tex
	pdflatex CgcssiRef
	bibtex CgcssiRef
	bibtex bu1
	bibtex bu2
	pdflatex CgcssiRef
	makeindex CgcssiRef.idx


refOld: CgcssiRef.tex
	latex CgcssiRef
	makeindex CgcssiRef.idx
	bibtex CgcssiRef
	latex CgcssiRef
	dvips -P pdf -G0 -o CgcssiRef.ps CgcssiRef.dvi
	ps2pdf CgcssiRef.ps
#	mv CgcssiRef.pdf $(WEBPAGE)/documentation

fig:
	latex fig.tex
	dvips -P pdf -G0 -o fig.ps fig.dvi
	ps2pdf fig.ps

clean: 
	-rm -f *.aux *.bbl *.blg *.ind *.ilg *.log *.idx *.toc
