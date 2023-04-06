# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'SurFEC'
copyright = '2022, Nicodemo Di Pasquale and Lorenzo Rovigatti'
author = 'Nicodemo Di Pasquale and Lorenzo Rovigatti'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinxarg.ext',
    'myst_parser',
    "sphinx_rtd_theme",
    "sphinxcontrib.bibtex"
]

bibtex_bibfiles = ['refs.bib']

intersphinx_mapping = {'python': ('https://docs.python.org/3/', None),
                       'numpy': ('https://numpy.org/doc/stable/', None),
                       'scipy': ('https://docs.scipy.org/doc/scipy/', None),
                       'matplotlib': ('https://matplotlib.org/stable/', None),
                       'sklearn' : ('https://scikit-learn.org/stable/', None)
                       }

# generate labels of heading anchors
myst_heading_anchors = 3

myst_enable_extensions = [
    "dollarmath",
    "amsmath"
    ]

napoleon_include_init_with_doc = True
autosummary_generate = True
autoclass_content = 'both'
autodoc_default_flags = ["show-inheritance", "members", "undoc-members"]
autodoc_default_options = {
    "show-inheritance": True, 
    "members": True, 
    "undoc-members": True
}

source_suffix = {
    '.rst' : 'restructuredtext',
    '.md' : 'markdown'
}

master_doc = 'index'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []
html_extra_path = []

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

def setup(app):
    pass
    
