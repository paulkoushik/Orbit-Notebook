.SUFFIXES: .pdf .png


.tex.pdf:
	pdflatex $<

.pdf.png:
	convert -density 300 $< $@

all: e7tree.png bfsdfs.png

clean:
	rm -f *.aux *.log

distclean: clean
	rm -f *.png


