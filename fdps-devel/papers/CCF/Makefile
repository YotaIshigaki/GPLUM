%.dvi : %.tex
	latex $*;
	latex $*;
	bibtex $*;
	latex $*;
%.pdf : %.dvi
	dvipdfmx  $*

SRCS= ccfpaper.tex
INSTALLDIR=/usr2/makino/src/fdps/fdps/papers/CCF
ccfpaper.pdf: $(SRCS)

install: $(SRCS)
	rsync -avprog $(SRCS) $(INSTALLDIR)
