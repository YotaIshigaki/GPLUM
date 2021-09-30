TRT=sc15_fdps
LTX=latex
BIB=bibtex
DVI=dvipdfmx

pdf:
	make tex
	make bib
	make tex
	make tex
	make dvi
tex:
	$(LTX) $(TRT).tex
bib:
	$(BIB) $(TRT).aux
dvi:
	$(DVI) $(TRT).dvi

clean:
	rm -f *~
distclean:
	make clean
	rm -f $(TRT).aux $(TRT).dvi $(TRT).log $(TRT).pdf $(TRT).blg $(TRT).bbl
