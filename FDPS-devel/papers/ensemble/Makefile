TARGET=ensemble-1
STY=ensemble.sty
EPS=$(shell ls *.eps)
TEX=platex
DVIPDF=dvipdfmx


.SUFFIXES: .tex .dvi

all:$(TARGET).dvi

$(TARGET).dvi: $(TARGET).tex $(EPS) $(STY)

.tex.dvi:
	$(TEX) $<
	$(TEX) $<
.pdf.dvi:
	$(DVIPDF) $<


tar:
	tar cvzf $(TARGET).tar.gz $(STY) $(TARGET).* $(EPS) makefile

pdf: $(TARGET).dvi
	$(DVIPDF) $(TARGET)

clean:
	rm -f $(TARGET).dvi $(TARGET).aux $(TARGET).log  $(TARGET).toc
	rm -f $(TARGET).pdf $(TARGET).tar.gz 
	rm -f *.pbm *.bmc
