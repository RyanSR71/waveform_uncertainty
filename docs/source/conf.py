# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'GWCorrect'
copyright = '2025, Ryan Johnson - No Rights Reserved'
author = 'Ryan Johnson'

release = 'beta'
version = '0.18.2'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'nbsphinx',
    'sphinx.ext.autosectionlabel',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'
pygments_style = "sphinx"

# -- Options for EPUB output
epub_show_urls = 'footnote'
