all:	aqm.tex arrayQualityMetrics.tex
	texi2dvi --pdf --clean arrayQualityMetrics.tex
	texi2dvi --pdf --clean aqm.tex 
	make clean

clean:
	rm -rf Report* Rplots.pdf aqm-* arrayQualityMetrics-* *.tex *.aux *.bbl *.blg *.log *.out *.toc 


