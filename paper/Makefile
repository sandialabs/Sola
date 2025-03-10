.PHONY: arxiv main toms clean help


arxiv:
	latexmk -pdf -quiet -g main.tex


main: arxiv


toms:
	latexmk -pdf -quiet -g main-toms.tex


clean:
	latexmk -c


help:
	@echo "make main  -> main.pdf      (for arXiv)"
	@echo "make arxiv -> main.pdf      (for arXiv)"
	@echo "make toms  -> main-toms.pdf (for journal)"
