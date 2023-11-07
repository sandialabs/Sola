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

```bash
python3 -m pip install -r requirements.txt
make html
```

To see the generated documentation, open `build/html/index.html` in a browser.

## Contributing

Documentation pages can be written in [reStructuredText](https://github.com/ralsina/rst-cheatsheet/blob/master/rst-cheatsheet.rst) (`.rst`) or (preferred) in [Markdown](https://myst-parser.readthedocs.io/en/latest/index.html) (`.md`).

### Useful MyST Links

- [MyST Markdown Guide](https://mystmd.org/guide/quickstart-myst-markdown)
- [Admonitions and Tabs](https://myst-parser.readthedocs.io/en/latest/syntax/admonitions.html) and [Admonitions in the Book Theme](https://sphinx-book-theme.readthedocs.io/en/stable/reference/kitchen-sink/admonitions.html)
- [Tables](https://myst-parser.readthedocs.io/en/latest/syntax/tables.html)
- [Math and Equations](https://myst-parser.readthedocs.io/en/latest/syntax/math.html)
- [Cross References](https://myst-parser.readthedocs.io/en/latest/syntax/cross-referencing.html)
- [Theorems, Proofs, and Algorithms](https://sphinx-proof.readthedocs.io/en/latest/syntax.html)

### Managing Dependencies

If more Python packages (e.g., sphinx extensions) are needed,

- add them (with version numbers) to [`requirements.in`](./requirements.in)
- install `pip-compile`: `python -m pip install pip-tools`
- run `make requirements` to update [`requirements.txt`](./requirements.txt).

**Note**: we were using the most up-to-date `sphinx==7.2.6`, but since it conflicts with `sphinx-book-theme==1.0.1`, we currently have `sphinx==6.2.1`.

### Documentation Organization

Documentation source files are contained in the [**source/**](./source/) directory, which has the following structure.

```python
source/
├── conf.py                         # Sphinx configuration file.
├── index.md                        # Root of documentation (home page).
├── _static/                        # CSS/HTML configuration.
├── problems/                       # Descriptions of outer-loop problems we treat.
│   ├── bayesinversion.md
│   ├── hdsa.md
│   ├── oed.md
│   └── optimization.md
├── sabl/                           # Documentation for SABL (MATLAB).
|   ├── about.md                    # SABL landing page.
|   ├── installation.md             # How to install SABL.
|   ├── anatomy.md                  # Table of contents for code structure guide.
|   ├── anatomy                     # Folder of code structure guides.
|   │   ├── bayesinversion.md
|   │   ├── hdsa.md
|   │   ├── oed.md
|   │   └── optimization.md
|   ├── examples.md                 # Table of contents for examples.
|   └── examples/                   # Folder of examples.
|       ├── example1.md
|       ├── example2.md
|       └── ...
├── mrhyde/                         # Documentation for MrHyDE (C++).
|   ├── about.md                    # MrHyDE landing page.
|   ├── installation.md             # How to install MrHyDE.
|   ├── primer.md                   # Transitioning from MATLAB to C++.
|   ├── anatomy.md                  # Table of contents for code structure guide.
|   ├── anatomy/                    # Folder of code structure guides.
|   │   ├── bayesinversion.md
|   │   ├── hdsa.md
|   │   ├── oed.md
|   │   └── optimization.md
|   ├── examples.md                 # Table of contents for examples.
|   └── examples/                   # Folder of examples.
|       ├── example1.md
|       ├── example2.md
|       └── ...
└── appendix/                       # Folder of appendices.
    ├── contributing.md             # How to contribute.
    ├── fundamentals.md             # Brief overview of mathematical prerequisites.
    └── notation.md                 # Index of notation.
```
