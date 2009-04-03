all: arrayQualityMetrics aqm clean

arrayQualityMetrics: arrayQualityMetrics.tex
	cp -p ../scripts/arrayQualityMetrics.pdf .

aqm: aqm.tex
	cp -p ../scripts/aqm.pdf .

clean:
	rm -f *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex Rplots.pdf


