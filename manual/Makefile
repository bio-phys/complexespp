MAIN=manual

all: $(MAIN).pdf

$(MAIN).pdf: src/*.tex  src/$(MAIN).sty
	cd src ; \
	latexmk -pdf -g $(MAIN).tex -shell-escape
	mv src/$(MAIN).pdf .


.PHONY: clean

clean:
	git clean -fdx
