# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'CorticalSim 3D'
copyright = '2026, Eva Deinum, Bandan Chakrabortty, Jaro Camphuijsen, Alex Hadjiivanov'
author = 'Eva Deinum, Bandan Chakrabortty, Jaro Camphuijsen, Alex Hadjiivanov'
release = ''

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['myst_parser', 'breathe', 'sphinx_rtd_theme']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Add breathe extension for c++ docs in sphinx with doxygen
breathe_projects = {"myproject": "../build/docs/doxygen/xml/"}
breathe_default_project = "myproject"
