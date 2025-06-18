# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
from urllib.request import urlopen
from pathlib import Path
from recommonmark.parser import CommonMarkParser

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'kinnex-docs-external'
copyright = '2024, Methods Dev Lab'
author = 'mdl@broadinstitute.org'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# -- General configuration ---------------------------------------------------

source_parsers = {
    '.md': CommonMarkParser,
}
source_suffix = ['.rst', '.md']

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["sphinx_design", "myst_nb", "sphinxemoji.sphinxemoji", "sphinx_thebe"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3.8", None),
    "mdl": ("https://github.com/MethodsDev/kinnex-documentation-external", None),
    "pst": ("https://pydata-sphinx-theme.readthedocs.io/en/latest/", None),
}
nitpick_ignore = [
    ("py:class", "docutils.nodes.document"),
    ("py:class", "docutils.parsers.rst.directives.body.Sidebar"),
]

suppress_warnings = ["myst.domains", "ref.ref"]

numfig = True

myst_enable_extensions = [
    "dollarmath",
    "amsmath",
    "deflist",
    # "html_admonition",
    "html_image",
    "colon_fence",
    # "smartquotes",
    # "replacements",
    # "linkify",
    # "substitution",
]

myst_url_schemes = ("http", "https", "mailto")

nb_execution_mode = "force"
nb_execution_allow_errors=True 

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"
html_logo = "_images/bcl_logo.png"
html_favicon = "_images/mdl_logo.png"

html_theme_options = {
    #'body_max_width' : 'none',
    'page_width': 'auto',
    "repository_url": "https://github.com/MethodsDev/kinnex-documentation-external",
    "repository_branch": "main",
    "launch_buttons": {
        "binderhub_url": "https://notebooks.gesis.org/binder/",
        "colab_url": "https://colab.research.google.com/",
        "deepnote_url": "https://deepnote.com/",
        "notebook_interface": "jupyterlab",
        "thebe": True,
    },
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/MethodsDev/kinnex-documentation-external",
            "icon": "fa-brands fa-github",
        }
    ],
    "use_edit_page_button": True,
    "use_source_button": True,
    "use_issues_button": True,
    "use_repository_button": True,
    "use_download_button": True,
    "use_sidenotes": True
}

thebe_config = {
   "codemirror-theme": "material-palenight"
}

def setup(app):
    app.add_css_file('my_theme.css')

html_static_path = ['_static']