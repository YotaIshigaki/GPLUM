TRT=sc15_memo
LTX=platex
DVI=dvipdfmx

pdf:
	make tex
	make dvi
tex:
	$(LTX) $(TRT).tex
	$(LTX) $(TRT).tex
dvi:
	$(DVI) $(TRT).dvi

clean:
	rm -f *~
distclean:
	make clean
	rm -f $(TRT).aux $(TRT).dvi $(TRT).log $(TRT).toc $(TRT).out
