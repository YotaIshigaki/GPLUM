#
# Makefile for acmart package
#
# This file is in public domain
#
# $Id: Makefile,v 1.10 2016/04/14 21:55:57 boris Exp $
#

PACKAGE=acmart



PDF = fdps-sunway.pdf

all:  ${PDF}

INSTALLDIR=/home/makino/src/fdps/fdps/papers/sc17

SRCS=fdps-sunway.tex fdps-sunway.bbl Makefile
UTILS= acmart.cls ACM-Reference-Format.bst

install: $(SRCS) $(UTILS) 
	rsync -avprog  $(SRCS)  $(UTILS)  $(INSTALLDIR) 

# %.pdf:  %.dtx   $(PACKAGE).cls
# 	pdflatex $<
# 	- bibtex $*
# 	pdflatex $<
# 	- makeindex -s gind.ist -o $*.ind $*.idx
# 	- makeindex -s gglo.ist -o $*.gls $*.glo
# 	pdflatex $<
# 	while ( grep -q '^LaTeX Warning: Label(s) may have changed' $*.log) \
# 	do pdflatex $<; done


# acmguide.pdf: $(PACKAGE).dtx $(PACKAGE).cls
# 	pdflatex -jobname acmguide $(PACKAGE).dtx
# 	- bibtex acmguide
# 	pdflatex -jobname acmguide $(PACKAGE).dtx
# 	while ( grep -q '^LaTeX Warning: Label(s) may have changed' acmguide.log) \
# 	do 	pdflatex -jobname acmguide $(PACKAGE).dtx; done

# %.cls:   %.ins %.dtx  
# 	pdflatex $<

%.pdf:  %.tex  # $(PACKAGE).cls ACM-Reference-Format.bst
	env TEXINPUTS=""  pdflatex $<
#	- bibtex $*
	env TEXINPUTS="" pdflatex $<
	env TEXINPUTS="" pdflatex $<
	while ( grep -q '^LaTeX Warning: Label(s) may have changed' $*.log) \
	do env TEXINPUTS="" pdflatex $<; done

fdps-sunway.pdf: fdps-sunway.tex

