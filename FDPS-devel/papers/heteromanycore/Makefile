PDF = hetero-fdps.pdf

all:  ${PDF}

OTHERFILES =  hetero.sty

FIGURES= 

%.dvi : %.tex
	platex $<; platex $<; bibtex $< ;platex $<; platex $<; 
%.pdf : %.dvi
	dvipdfmx  $<

hetero-fdps.dvi: hetero-fdps.tex $(FIGURES) $(BIBFILES) $(OTHERFILES)

