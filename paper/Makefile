.PHONY: all arxiv toms clean help


all: arxiv toms


arxiv:
	latexmk -pdf -quiet -g --shell-escape main-arxiv.tex


toms:
	latexmk -pdf -quiet -g --shell-escape main-toms.tex


clean:
	latexmk -c


help:
	@echo "make arxiv -> main-arxiv.pdf  (for arXiv)"
	@echo "make toms  -> main-toms.pdf   (for journal)"
	@echo "make clean -> clean extra latex files"
