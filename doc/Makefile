LATEXMK=latexmk -pdf

%.pdf : %.tex
	$(LATEXMK) $<

.PHONY: clean

# remove intermediate files
clean:
	@rm -f *.aux
	@rm -f *.fdb_latexmk
	@rm -f *.fls
	@rm -f *.log
	@rm -f *.nav
	@rm -f *.out
	@rm -f *.snm
	@rm -f *.synctex.gz
	@rm -f *.toc
	@rm -f *.dvi
	@rm -f *.bbl
	@rm -f *.bcf
	@rm -f *.blg
	@rm -f *.run.xml
