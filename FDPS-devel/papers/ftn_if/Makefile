# Target objects
DOC = ms
# Shell commands
LATEX = platex -shell-escape
DVIPDFMX = dvipdfmx

all:
	# Make *.dvi
	$(LATEX) $(DOC) 
	$(LATEX) $(DOC) 
	$(LATEX) $(DOC) 
	# Make *.pdf
	$(DVIPDFMX) $(DOC)
	# Open [Mac OSX only]
	#open $(DOC).pdf

clean:
	rm -f *.aux *.dvi *.log *.toc *.out

distclean:
	make clean
	rm -f *.pdf
