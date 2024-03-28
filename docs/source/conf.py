# conf.py
"""Configuration file for the Sphinx documentation builder.

This file only contains a selection of the most common options.
For a full list see the documentation:
https://www.sphinx-doc.org/en/master/usage/configuration.html.
"""

import os
import time


# -- Project information -----------------------------------------------------

project = "WOLF"
copyright = f"{time.strftime('%Y')} Sandia National Laboratories"
author = "Joseph Hart, Shane A. McQuarrie, and Bart van Bloemen Waanders"

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
    # "sphinx_tippy",  # Previews when hovering over links.
    "sphinx_togglebutton",  # Dropdowns.
    "sphinxcontrib.matlab",  # MATLAB API generation.
    # "sphinxcontrib.bibtex",  # LaTeX citations.
    "sphinx.ext.autodoc",  # API generation.
    "sphinx.ext.napoleon",  # Better docstring formatting.
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
html_theme_options = {
    "path_to_docs": "docs/",
    "home_page_in_toc": True,
    "repository_url": "https://cee-gitlab.sandia.gov/joshart/sabl",
    # "use_issues_button": True,    # Not supported for GitLab
    "use_download_button": False,
    "use_fullscreen_button": True,
    "use_repository_button": True,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = [
    "admonish.css",  # Custom admonition :class: abstract
    "colors.css",  # Explicit text coloring.
    "rightalign.css",  # Right align equation labels.
    # "tippy.css",  # Correct theme for link previews.
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
macros = {
    "trp": r"{^{\mathsf{T}}}",
    "invtrp": r"{^{-\mathsf{T}}}",
    "grad": [r"\nabla_{\!#1}\,", 1],
    "ddt": r"\frac{\textrm{d}}{\textup{d}t}",
    "J": r"\mathcal{J}",
    "N": r"\mathbb{N}",
    "R": r"\mathbb{R}",
}
macros.update({x: rf"\mathbf{{{x}}}" for x in "bcfghquvwxyzABCDHIQSVWXYZ01"})
macros.update(
    {
        f"bf{x}": rf"\boldsymbol{{\{x}}}"
        for x in ["lambda", "mu", "gamma", "Gamma", "Phi", "Sigma", "Psi"]
    }
)
mathjax3_config = {"tex": {"macros": macros}}

# LaTeX citations.
# bibtex_bibfiles = ["references.bib"]

# MATLAB API generation.
autodoc_member_order = "bysource"
matlab_auto_link = "basic"
matlab_class_signature = True
matlab_short_links = True
matlab_show_property_default_value = True
matlab_src_dir = os.path.abspath(os.path.join("..", "..", "src"))
primary_domain = "mat"
