PROJECT = '2021_SeaIceDeformation'
MAIN = SeaIceDeformation_Methods
TEX_SOURCES = Makefile \
              $(MAIN).tex

FIGURES := $(shell find figures/* -type f)

SHELL=/bin/bash
DATE = $(shell date +"%d%b%Y")
OPT = --interaction=nonstopmode

all: $(MAIN).pdf

final: final.pdf
	$(MAKE) clean

$(MAIN).pdf: $(TEX_SOURCES) $(FIGURES)
	pdflatex $(MAIN)
	bibtex $(MAIN)
	pdflatex $(MAIN)
	pdflatex $(MAIN)
	@#latexmk -pdf -pvc -pdflatex="pdflatex $(OPT)" $(MAIN) 

once: 
	pdflatex $(MAIN)

clean: 
	-rm -f $(MAIN).{pdf,log,blg,bbl,aux,out,toc,idx,bcf,run.xml,ind,ilg,fls,fdb_latexmk} $(MAIN).pdf


.PHONY: clean all
