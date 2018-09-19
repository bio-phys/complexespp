MAIN=manual

all: $(MAIN).pdf

$(MAIN).pdf: *.tex  $(MAIN).sty
	latexmk -pdf -g $(MAIN).tex -shell-escape


.PHONY: clean

clean:
	latexmk -c
	rm -f $(MAIN).acn $(MAIN).acf $(MAIN).alg $(MAIN).bbl $(MAIN).glg
	rm -f $(MAIN).glo $(MAIN).gls $(MAIN).ist $(MAIN).acr
	rm -rf _minted-$(MAIN)
