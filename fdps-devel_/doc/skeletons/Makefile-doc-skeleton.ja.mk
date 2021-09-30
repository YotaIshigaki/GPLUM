#
# Makefile-doc-skelton.ja.mk
#
# Skelton Makefile for FDPS Japanese document
# Shell commands
LATEX = platex -shell-escape
DVIPDFMX = dvipdfmx 

# Pattern rules
%.pdf : %.dvi
	$(DVIPDFMX)  $<

# Target objects
DOC = $(basename $(wildcard doc*.tex))
TEX = $(DOC).tex
DVI = $(DOC).dvi
PDF = $(DOC).pdf

# Source files
SRC = $(wildcard *.tex)

all: $(PDF)
	# Copy 
	cp $(DOC).pdf ../../../
	# Open [Mac OSX only]
	#open $(DOC).pdf

$(DVI): $(SRC)
	$(LATEX)  $(DOC); $(LATEX) $(DOC); $(LATEX) $(DOC);

clean:
	rm -f *.aux *.dvi *.log *.toc *.out

distclean: clean
	rm -f $(PDF)
