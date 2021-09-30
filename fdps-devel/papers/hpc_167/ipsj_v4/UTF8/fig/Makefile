# Shell commands
LATEX = platex -shell-escape
PDF2EPS = pdf2eps
DVIPDFMX = dvipdfmx

# Pattern rules
%.eps : %.pdf
	$(PDF2EPS) 1 $(basename $<)

%.pdf : %.dvi
	$(DVIPDFMX) $(basename $<)

%.dvi : %.tex
	$(LATEX) $(basename $<)

# Target objects
SRC = $(wildcard *.tex)
DVI = $(SRC:tex=dvi)
PDF = $(SRC:tex=pdf)
EPS = $(SRC:tex=eps)

all: $(EPS)

clean:
	rm -f *.aux *.dvi *.log *.toc *.out

distclean:
	make clean
