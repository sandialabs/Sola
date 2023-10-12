# SABL Documentation

This folder contains documentation for SABL and its larger mission of connecting users to high-performance PDE-constrained optimization and other outer-loop tasks.

## Contents

- [`main.tex`](./main.text): technical summary document, not part of the web-based documentation.
- [`Makefile`](./Makefile): Makefile for the documentation. Run `make help` for options.
- [**source/**](./source/): source files for the web-based documentation.
- **build/**: directory containing compiled HTML files. This folder is not tracked by git.


## Building the documentation

**Summary document**

```
$ make main.pdf
```

**Web-based Documentation**

```
$ make html
```

To see the generated documenation, open **TODO** in a browser.

