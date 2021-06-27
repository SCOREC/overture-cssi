
all: ogenDoc
doc: ogenDoc

OGEN = $(OvertureCheckout)/ogen

ogenUpdate.tex: $(OGEN)/checkOverlap.C computeOverlap.C
	extract.p $(OGEN)/checkOverlap.C $(OGEN)/computeOverlap.C
changeParametersInclude.tex: $(OGEN)/changeParameters.C
	extract.p $(OGEN)/changeParameters.C
ogenUpdateInclude.tex: $(OGEN)/checkOverlap.C $(OGEN)/computeOverlap.C
	extract.p $(OGEN)/checkOverlap.C $(OGEN)/computeOverlap.C


docIncludeFiles = changeParametersInclude.tex ogenInclude.tex ogenUpdateInclude.tex

ogenDoc: $(docIncludeFiles)
	pdflatex ogen
	bibtex ogen
	makeindex ogen.idx

ogenDocOld: $(docIncludeFiles)
	latex ogen
	makeindex ogen.idx
	latex ogen
	dvips -P pdf -G0 -o ogen.ps ogen.dvi
	ps2pdf ogen.ps
	cp ogen.pdf $(WEBPAGE)/documentation

ogenhtml: $(docIncludeFiles)
	lh.p ogen.tex
	latex ogenHTML
	bibtex ogenHTML
	makeindex ogenHTML.idx
	latex ogenHTML
	latex2html -local_icons ogenHTML
	rm -rf $(WEBPAGE)/documentation/ogenHTML
	mkdir $(WEBPAGE)/documentation/ogenHTML
	cp ogenHTML/*.* $(WEBPAGE)/documentation/ogenHTML
	rm -rf ogenHTML


cartoon: 
	pdflatex cartoon


