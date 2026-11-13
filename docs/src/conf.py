# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import tomllib
from datetime import datetime

# -- Project information -----------------------------------------------------

# The version is read from Cargo.toml so the documentation always matches the crate.
cargo_toml = os.path.join(os.path.dirname(__file__), '..', '..', 'Cargo.toml')
with open(cargo_toml, 'rb') as f:
    cargo = tomllib.load(f)['package']

project = 'gTrim'
author = ', '.join(cargo['authors'])
copyright = f'{datetime.now().year}, {author}'

release = cargo['version']
version = cargo['version']

github_url = 'https://github.com/cmc-aau/gTrim'

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# -- General configuration ---------------------------------------------------

extensions = ['myst_parser', 'sphinx_sitemap']

templates_path = ['_templates']
exclude_patterns = []

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

html_extra_path = ['robots.txt']

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    'logo_only': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': True,
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
}

html_static_path = ['_static']

html_context = {
    'display_github': True,
    'github_user': 'cmc-aau',
    'github_repo': 'gTrim',
    'github_version': 'main/docs/src/',
}

html_css_files = ['theme_overrides.css']
html_js_files = ['custom.js']

html_logo = '_static/gtrim_stacked.svg'
html_favicon = '_static/gtrim_favicon.png'

# Sitemap settings
html_baseurl = 'https://cmc-aau.github.io/gTrim/'
sitemap_url_scheme = '{link}'


# Replace the version placeholder in all pages, including code blocks, with the version from Cargo.toml.
VERSION_PLACEHOLDER = '{{ gtrim_version }}'


def replace_version_placeholder(app, docname, source):
    source[0] = source[0].replace(VERSION_PLACEHOLDER, release)


def setup(app):
    app.connect('source-read', replace_version_placeholder)
