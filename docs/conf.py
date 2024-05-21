# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

import os
import sys
import re
import sphinx
import mock

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

# project root
sys.path.insert(0, os.path.abspath('../'))
#sys.path.append('/home/jramos/PycharmProjects/chiripa') 
#sys.path.append('/opt/pycharm-2020.1/plugins/python/helpers/rest_runners')
#sys.path.append('/home/jramos/PycharmProjects/chiripa') 
#sys.path.append('/opt/pycharm-2020.1/plugins/python/helpers/pycharm_display')
#sys.path.append('/opt/pycharm-2020.1/plugins/python/helpers')
#sys.path.append('/usr/lib/python38.zip')
#sys.path.append('/usr/lib/python3.8')
#sys.path.append('/usr/lib/python3.8/lib-dynload')
#sys.path.append('/home/jramos/PycharmProjects/sandbox_chiripa/lib/python3.8/site-packages')
#sys.path.append('/home/jramos/PycharmProjects/sandbox_chiripa/lib/python3.8/site-packages/PolyAnaGro-0.1-py3.8-linux-x86_64.egg')
#sys.path.append('/home/jramos/PycharmProjects/sandbox_chiripa/lib/python3.8/site-packages/Chiripa-0.1-py3.8-linux-x86_64.egg')
#sys.path.append('/opt/pycharm-2020.1/plugins/python/helpers/pycharm_matplotlib_backend')
#print(sys.path)

# -- Project information -----------------------------------------------------
project = 'Chiripa -- CHI inteRactIon PArameter.'
copyright = '2019, J.Ramos'
author = 'J.Ramos'
# The full version, including alpha/beta/rc tags
version = '0.1'
release = '0.1'

# -- General configuration ---------------------------------------------------
extensions = [
        'sphinx.ext.autodoc',
        'sphinx.ext.intersphinx',
        'sphinx.ext.autosummary',
        'sphinx.ext.napoleon',
        'sphinx.ext.mathjax',
        'sphinx.ext.viewcode'
]

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = False
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = False
napoleon_type_aliases = None

#autosummary_generate = True

intersphinx_mapping = {'python': ('https://docs.python.org/3', None),
                       'numpy': ('https://docs.scipy.org/doc/numpy', None)}

autodoc_docstring_signature = True

sphinx_ver = tuple(map(int, sphinx.__version__.split('.')))
if sphinx_ver < (1,8,0):
    autodoc_default_flags = ['inherited-members']
else:
    autodoc_default_options = {'inherited-members': None}

#autodoc_mock_imports = ['numpy','h5py', 'networkx', 'scipy', 'matplotlib']
MOCK_MODULES = ['pandas', 'networkx', 'matplotlib', 'numpy', 'scipy', 'h5py', 'ext_libc', 'ext_libc.c_discriminant']
for mod_name in MOCK_MODULES:
    sys.modules[mod_name] = mock.Mock()

templates_path = ['_templates']

source_suffix = '.rst'

master_doc = 'index'

todo_include_todos = True

language = None

# Relative to source directory
exclude_patterns = ['_build', '_templates']

pygments_style = 'sphinx'

# -- HTML Theme ---------------------------------------------------
html_theme = 'sphinx_rtd_theme'
#html_theme = 'pydata_sphinx_theme'
#html_theme = 'bizstyle'
#html_theme = 'pyramid'
#html_theme = 'nature'
#html_theme = 'haiku'

#html_static_path = ['_static']

# This is required for the alabaster theme
# refs: http://alabaster.readthedocs.io/en/latest/installation.html#sidebars
html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html',  # needs 'show_related': True theme option to display
        'searchbox.html',
        'donate.html',
    ]
}


htmlhelp_basename = 'chiripa-doc'



# -- Options for LaTeX output ---------------------------------------------
latex_elements = {
# The paper size ('letterpaper' or 'a4paper').
#'papersize': 'letterpaper',

# The font size ('10pt', '11pt' or '12pt').
#'pointsize': '10pt',

# Additional stuff for the LaTeX preamble.
#'preamble': '',

# Latex figure (float) alignment
#'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'chiripa.tex', 'CHIRIPA Documentation',
     'Javier Ramos, IEM-CSIC', 'manual'),
]

# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'chiripa', 'CHIRIPA Documentation',
     [author], 1)
]

# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'CHIRIPA', 'CHIRIPA Documentation',
     author, 'CHIRIPA', 'One line description of project.',
     'Miscellaneous'),
]

ipython_mplbackend = None;
