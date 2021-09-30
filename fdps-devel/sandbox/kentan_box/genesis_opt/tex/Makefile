FILE=main
PLATEX=uplatex
TARGET=$(FILE).pdf

all	:$(TARGET)

$(TARGET)	:$(FILE).dvi 
	dvipdfmx $(FILE).dvi

$(FILE).dvi	:$(FILE).tex
	$(PLATEX) $(FILE).tex
	bibtex 	$(FILE)
	$(PLATEX) $(FILE).tex

open: all
	open $(TARGET)

back	:
	mkdir bkup
	cp -r mypaper.tex Makefile bkup/.
	tar zcvf bkup$(shell date +"%Y%m%d").tar.gz bkup
	rm -r bkup

clean	:
	rm $(FILE).dvi $(FILE).pdf $(FILE).aux $(FILE).bbl
