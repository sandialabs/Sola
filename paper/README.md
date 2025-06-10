# Paper

Draft for a software paper on the MATLAB/SABL to C++/MrHyDE framework, featuring HdsaLib.

Potential target: [TOMS](https://dl.acm.org/journal/toms), in the **Research Paper** category.

## Paper compilation

The paper uses the `minted` package to render nice code blocks.
This requires the python package `Pygmentize` (`pip install Pygmentize`) and to use the `-shell-escape` flag when compiling with `pdflatex` or `latexmk`.
If you have this set up, use this import line in `config.tex`.

```latex
\usepackage[finalizecache,cachedir=_minted]{minted}
```

Otherwise, use this line and **do not add any new `minted` environments to the tex.**

```latex
\usepackage[frozencache,cachedir=_minted]{minted}
```
