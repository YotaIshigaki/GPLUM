# Target objects
SRC_DISPATCH_INDEX = src_dispatch_index
SRC_DISPATCH = src_dispatch
SRC_RETRIEVE = src_retrieve
SRC_REUSE = src_reuse
# Shell commands
LATEX = latex -shell-escape
DVIPS = dvips
PS2EPS = ps2epsi

all:
	make src_dispatch
	make src_dispatch_index
	make src_retrieve
	make src_reuse
	make distclean

src_dispatch:
	$(LATEX)  $(SRC_DISPATCH).tex
	$(DVIPS)  $(SRC_DISPATCH).dvi
	$(PS2EPS) $(SRC_DISPATCH).ps $(SRC_DISPATCH).eps

src_dispatch_index:
	$(LATEX)  $(SRC_DISPATCH_INDEX).tex
	$(DVIPS)  $(SRC_DISPATCH_INDEX).dvi
	$(PS2EPS) $(SRC_DISPATCH_INDEX).ps $(SRC_DISPATCH_INDEX).eps

src_retrieve:
	$(LATEX)  $(SRC_RETRIEVE).tex
	$(DVIPS)  $(SRC_RETRIEVE).dvi
	$(PS2EPS) $(SRC_RETRIEVE).ps $(SRC_RETRIEVE).eps

src_reuse:
	$(LATEX)  $(SRC_REUSE).tex
	$(DVIPS)  $(SRC_REUSE).dvi
	$(PS2EPS) $(SRC_REUSE).ps $(SRC_REUSE).eps

clean:
	rm -f *.aux *.dvi *.log *.toc *.out

distclean:
	make clean
	rm -f *.dvi
