#!/bin/sh

for FILE in domain_c*.eps
            
do
    IN=${FILE%.eps}
    ps2pdf -dEPSCrop $IN.eps $IN.pdf
    pdftops -eps $IN.pdf $IN.eps
    rm $IN.pdf
done
