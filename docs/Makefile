THEORY=theory
ARCHITECTURE=architecture

all: $(THEORY).pdf $(ARCHITECTURE).pdf

$(THEORY).pdf: theory/*.tex  theory/$(THEORY).sty
	cd theory ; \
	latexmk -pdf -g $(THEORY).tex -shell-escape
	mv theory/$(THEORY).pdf .


$(ARCHITECTURE).pdf: architecture/*.tex  architecture/$(ARCHITECTURE).sty
	cd architecture ; \
	latexmk -pdf -g $(ARCHITECTURE).tex -shell-escape
	mv architecture/$(ARCHITECTURE).pdf .


.PHONY: clean

clean:
	git clean -fdx
