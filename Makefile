##################################################################
# Makefile for LaTeX
##################################################################
# Use:
# make
# make clean
TEX:=$(shell ls *.tex)
FILES = *.tex Makefile *.bst *.pdf *.bib
FOLDER = ./home/saul/sauld@cimat.mx/UNISON/Articles/FSDEhiv-aidsModel/git_fbmSDEModels/versions
TEXLIST = TexList.txt
FIGLIST = FigList.txt
CHORILIST = ChoriList.txt
OTHER = *~ *.aux *.dvi *.toc *.bbl *.blg *.out *.thm *.ps *.idx *.ilg *.ind *.tdo *.pdf *.tar.gz *.log *.spl *.synctex.gz
NAMETAR1:= $(shell date '+%Y%b%d_%H_%M')
NAMETAR = "$(NAMETAR1)fbm_sde_hiv_model.tar.gz"
NAMETARTEX = "$(NAMETAR1)fbm_sde_hiv_model.tar.gz"
NAMETARFIG = "$(NAMETAR1)fbm_sde_hiv_model.tar.gz"
NAMETARCHORI = "$(NAMETAR1)fbm_sde_hiv_model.tar.gz"
pdflatex: fbm_sde_model.tex
	pdflatex --synctex=1 fbm_sde_model.tex
	bibtex fbm_sde_model.aux
	pdflatex --synctex=1 fbm_sde_model.tex
	pdflatex --synctex=1 fbm_sde_model.tex

clean:
	rm -f $(OTHER)
#
tar:
	tar -cvf $(NAMETAR) -T $(TEXLIST)
#
tarFig:
	(cd ..&& tar -cvf $(NAMETARFIG) -T $(FIGLIST))
	
#
tarTex:
	(cd ..&& tar -cvf $(NAMETARTEX) -T $(TEXLIST))

tarChori:
	tar -cvf $(NAMETARCHORI) -T $(CHORILIST)
	
Chori:	Chorizo.sh
	sh Chorizo.sh
	pdflatex --synctex=1 fbm_sde_model.tex
	bibtex fbm_sde_model.aux
	pdflatex --synctex=1 fbm_sde_model.tex
	pdflatex --synctex=1 fbm_sde_model.tex
