# SABL Documentation

This folder contains documentation for `SABL` and its larger mission of connecting users to high-performance PDE-constrained optimization and other outer-loop tasks.

## Contents

- [`Makefile`](./Makefile): Makefile for the documentation. Run `make help` for options.
- [**overview/**](./overview/): LaTex project overview document.
- [**source/**](./source/): source files for the web-based documentation.
- **build/**: directory containing compiled HTML files. This folder is not tracked by git.
- [`requirements.in`](./requirements.in): high-level requirements for compiling the web-based documentation.
- [`requirements.txt`](./requirements.txt): exhaustive requirements for compiling the web-based documentation.


## Building the documentation

The documentation is built with the `sphinx` Python infrastructure, see [the sphinx documentation](https://docs.readthedocs.io/en/stable/intro/getting-started-with-sphinx.html).

```
$ python3 -m pip install -r .requirements.txt
$ make html
```

To see the generated documentation, open `build/html/index.html` in a browser.

## Contributing

Documentation pages can be written in [reStructuredText](https://github.com/ralsina/rst-cheatsheet/blob/master/rst-cheatsheet.rst) (`.rst`) or in [Markdown](https://myst-parser.readthedocs.io/en/latest/index.html) (`.md`).

### Managing Dependencies

If more Python packages are needed,
- add them (with version numbers) to [`requirements.in`](./requirements.in)
- install `pip-compile`: `python -m pip install pip-tools`
- run `make requirements` to update [`requirements.txt`](./requirements.txt).
