# conf.py
"""
Configuration file for the Sphinx documentation builder.

This file only contains a selection of the most common options.
For a full list see the documentation:
https://www.sphinx-doc.org/en/master/usage/configuration.html.
"""

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

import time


project = "WOLF"
copyright = f"{time.strftime('%Y')} Sandia National Laboratories"
author = "Joseph Hart, Shane McQuarrie, and Bart van Bloemen Waanders"

# The full version, including alpha/beta/rc tags
release = "0.0.1"

html_title = "WOLF 🐺"
html_short_title = "WOLF"
# html_logo = None
html_favicon = "../img/favicon.png"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings.
extensions = [
    "myst_parser",  # Write content in Markdown.
    "sphinx_copybutton",  # Copy code blocks.
    "sphinx_proof",  # Theorems, algorithms, etc.
    "sphinx_tippy",  # Previews when hovering over links.
    "sphinx_togglebutton",  # Dropdowns.
]

# Add any paths that contain templates here, relative to this directory.
# templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.
html_theme = "sphinx_book_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = [
    "colors.css",  # Explicit text coloring.
    "rightalign.css",  # Right align equation labels.
    "tippy.css",  # Correct theme for link previews.
]

# Explicit text coloring (ignoring theme).
rst_prolog = """
.. include:: <s5defs.txt>

"""

# -- Extensions --------------------------------------------------------------

myst_enable_extensions = [
    "amsmath",  # parse amsmath directly (align, etc.)
    "attrs_block",  # label and reference paragraphs, etc.
    "colon_fence",  # parse ::: delimiters.
    "dollarmath",  # Parse $ and $$ encapsulated math.
    "replacements",  # Convert (c) to ©, etc.
]
myst_dmath_double_inline = True  # $$ OK if no newline before and after.
suppress_warnings = ["myst.domains"]  # Suppress warning for using sphinx_proof

# LaTeX macros.
mathjax3_config = {
    "tex": {
        "macros": {
            "trp": r"{^{\mathsf{T}}}",
            "RR": r"\mathbb{R}",
        }
    }
}
