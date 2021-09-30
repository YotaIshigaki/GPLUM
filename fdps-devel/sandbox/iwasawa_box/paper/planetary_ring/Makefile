TRT=ms
LTX=platex
DVI=dvipdfmx
VIEW=evince
#VIEW=acroread
#VIEW=/cygdrive/c/Program\ Files\ \(x86\)/Adobe/Acrobat\ Reader\ DC/Reader/AcroRd32.exe

all:
	make tex
	make dvi
	make view

tex:
	$(LTX) $(TRT).tex
	$(LTX) $(TRT).tex

dvi:
	$(DVI) $(TRT).dvi

view:
	$(VIEW) $(TRT).pdf

clean:
	rm -f *~

distclean:
	make clean
	rm -f $(TRT).aux $(TRT).dvi $(TRT).log $(TRT).toc

