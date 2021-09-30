TARGET=plasma_fusion
STY=ensemble.sty
#EPS=$(shell ls *.eps)
EPS=  # JM made this change on Sep 4 2019
TEX=platex
DVIPDF=dvipdfmx

SRCS = $(TARGET).tex  intro.tex fdps_overview.tex

all: $(TARGET).pdf

$(TARGET).dvi: $(SRCS)  $(EPS) $(STY)


%.dvi : %.tex
	$(TEX) $<; $(TEX) $<;
%.pdf : %.dvi
	$(DVIPDF)  $<

tar:
	tar cvzf $(TARGET).tar.gz $(STY) $(TARGET).* $(EPS) makefile

#pdf: $(TARGET).dvi
#	$(DVIPDF) $(TARGET)

clean:
	rm -f $(TARGET).dvi $(TARGET).aux $(TARGET).log  $(TARGET).toc
	rm -f $(TARGET).pdf $(TARGET).tar.gz 
	rm -f *.pbm *.bmc
