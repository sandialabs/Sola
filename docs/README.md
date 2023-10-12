# SABL Documentation

This folder contains documentation for SABL and its larger mission of connecting users to high-performance PDE-constrained optimization and other outer-loop tasks.

## Contents

- [`main.tex`](./main.text): technical summary document, not part of the web-based documentation.
- [`Makefile`](./Makefile): Makefile for the documentation. Run `make help` for options.
- [**source/**](./source/): source files for the web-based documentation.
- **build/**: directory containing compiled HTML files. This folder is not tracked by git.
- [`.requirements.in`](./requirements.in): high-level requirements for compiling the web-based documentation.
- [`.requirements.txt`](./requirements.txt): exhaustive requirements for compiling the web-based documentation.


## Building the documentation

### Summary Document

```
$ make main.pdf
```

This method requires `latexmk`.

### Web-based Documentation

The documentation is built with the `sphinx` Python infrastructure, see [the sphinx documentation](https://docs.readthedocs.io/en/stable/intro/getting-started-with-sphinx.html).

```
$ python3 -m pip install -r .requirements.txt
$ make html
```

To see the generated documentation, open `build/html/index.html` in a browser.

#### Maintenance

If more Python packages are needed,
- add them (with version numbers) to `.requirements.in`
- install `pip-compile`: `python -m pip install pip-tools`
- run `make requirements`.
