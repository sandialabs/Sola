.PHONY: main


help:
	@echo "arxiv -> main-arxiv.pdf"
	@echo "toms  -> main-toms.pdf"


arxiv:
	latexmk -pdf -shell-escape main-arxiv.tex


toms:
	pip install --upgrade Pygments
	latexmk -pdf -quiet -Werror -g -shell-escape main-toms.tex
