"""Sphinx configuration for the dandi_acq4 documentation.

Build with:  uv run sphinx-build -b html docs docs/_build/html
"""

import sys
from pathlib import Path

# the package is "src" in the project directory
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

project = "dandi_acq4"
author = "Paul B. Manis"
copyright = "2022-2026, Paul B. Manis"
release = "0.3.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",  # numpy-style docstrings
    "sphinx.ext.viewcode",
    "myst_parser",  # Markdown pages
]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
myst_enable_extensions = ["dollarmath"]
exclude_patterns = ["_build"]

napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_use_ivar = True  # dataclass fields are also documented by autodoc
autodoc_member_order = "bysource"
autodoc_default_options = {"members": True, "undoc-members": True, "show-inheritance": True}
# display_acq4 needs the ephys package, which is not installed by default
autodoc_mock_imports = ["ephys"]

html_theme = "alabaster"
html_theme_options = {"page_width": "1100px", "sidebar_width": "240px"}
