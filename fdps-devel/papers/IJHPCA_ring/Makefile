TRT=ms
LTX=latex
BIB=bibtex
DVI=dvipdfmx
VIEW=evince
#VIEW=acroread
#VIEW=/cygdrive/c/Program\ Files\ \(x86\)/Adobe/Acrobat\ Reader\ DC/Reader/AcroRd32.exe

all:
	make tex
	make dvi
	#make copy
	#make view

tex:
	$(LTX) $(TRT).tex
	$(LTX) $(TRT).tex
	$(BIB) $(TRT)

dvi:
	$(DVI) $(TRT).dvi

copy:
	cp $(TRT).pdf ~/windows_dir/

view:
	$(VIEW) $(TRT).pdf

clean:
	rm -f *~

distclean:
	make clean
	rm -f $(TRT).aux $(TRT).dvi $(TRT).log $(TRT).toc

