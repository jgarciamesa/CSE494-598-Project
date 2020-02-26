LATEXMK=latexmk -pdf

%.pdf : %.tex
	$(LATEXMK) $<

.PHONY: clean

# remove intermediate files
clean:
	@rm -f *.{aux,fdb_latexmk,fls,log,nav,out,snm,synctex.gz,toc,dvi}
