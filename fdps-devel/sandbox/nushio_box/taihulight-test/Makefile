all: learn-sunway.pdf

learn-sunway.pdf: learn-sunway.tex */*/remark.tex
	./summarize-experiment.py
	platex --shell-escape learn-sunway.tex
	platex --shell-escape learn-sunway.tex
	dvipdfmx learn-sunway.dvi
