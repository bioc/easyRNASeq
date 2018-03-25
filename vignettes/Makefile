# Makefile to use knitr for package vignettes

## The Targets
## -----------
.PHONY: usage

all: pdf html clean

pdf:
	$(R_HOME)/bin/R CMD Sweave easyRNASeq.Rnw
	$(R_HOME)/bin/Rscript -e "tools::texi2pdf('easyRNASeq.tex')"

html:
	$(R_HOME)/bin/Rscript -e "rmarkdown::render('simpleRNASeq.Rmd')"

clean:
	@rm -rf *.tex *.bbl *.blg *.aux *.out *.log *.spl *.md

usage:
	@if [ -z $(R_HOME) ]; then \
	  echo "usage:\n \
		$(R_HOME)/bin/R CMD make \n \
		SWEAVE_OPTIONS=\"eval=FALSE\" $(R_HOME)/bin/R CMD make" ; \
	  exit 1 ; \
	fi
