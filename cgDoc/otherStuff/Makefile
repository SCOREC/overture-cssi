all = otherStuff
all: $(all);

MAPPING = $(OvertureCheckout)/mapping
GF = $(OvertureCheckout)/gf
OTHERSTUFF = $(OvertureCheckout)/otherStuff


otherStuffFiles = $(OTHERSTUFF)/OGgetIndex.C \
                  $(OTHERSTUFF)/getFromADataBase.C \
                  $(OTHERSTUFF)/OGPolyFunction.bC \
                  $(OTHERSTUFF)/OGFunction.C \
                  $(OTHERSTUFF)/OGTrigFunction.C \
                  $(OTHERSTUFF)/display.C \
                  $(OTHERSTUFF)/displayMask.C \
                  $(OTHERSTUFF)/Integrate.C \
                  $(OTHERSTUFF)/TridiagonalSolver.C \
                  $(MAPPING)/sPrintF.C \
                  $(OTHERSTUFF)/FortranIO.C \
                  $(OTHERSTUFF)/OGPulseFunction.C \
                  $(GF)/OvertureInit.C

otherStuff:
	extract.p $(otherStuffFiles)
	pdflatex otherStuff
	makeindex otherStuff


otherStuffOld:
	extract.p $(otherStuffFiles)
	latex otherStuff
	makeindex otherStuff
	latex otherStuff
	dvips -P pdf -G0 -o otherStuff.ps otherStuff.dvi 
	ps2pdf otherStuff.ps
	cp otherStuff.pdf $(WEBPAGE)/documentation

otherStuffhtml: 
	lh.p otherStuff.tex
	latex otherStuffHTML
	# bibtex otherStuffHTML
	makeindex otherStuffHTML.idx
	latex otherStuffHTML
	latex2html -local_icons otherStuffHTML
	rm -rf $(WEBPAGE)/documentation/otherStuffHTML
	mkdir $(WEBPAGE)/documentation/otherStuffHTML
	cp otherStuffHTML/*.* $(WEBPAGE)/documentation/otherStuffHTML
	rm -rf otherStuffHTML

cutCellsDoc:
	extract.p CutCells.C
	latex cutCells
	makeindex cutCells



