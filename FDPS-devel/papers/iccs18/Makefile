TRT=ms
LTX=latex
BIB=bibtex
DVI=dvipdfmx
#VIEW=acroread
VIEW=/cygdrive/c/Program\ Files\ \(x86\)/Adobe/Acrobat\ Reader\ DC/Reader/AcroRd32.exe

all:
	make tex
	make bib
	make tex
	make dvi

tex:
	$(LTX) $(TRT).tex
	$(LTX) $(TRT).tex

bib:
	$(BIB) $(TRT)
dvi:
	$(DVI) $(TRT).dvi

view:
	$(VIEW) $(TRT).pdf

clean:
	rm -f *~

distclean:
	make clean
	rm -f $(TRT).aux $(TRT).dvi $(TRT).log $(TRT).toc

