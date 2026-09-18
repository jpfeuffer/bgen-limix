"""scikit-build-core dynamic-metadata plugin.

Adds the nanobind-backend runtime dependency only when building in the
default "shared" nanobind split-mode (used for PyPI wheels). The conda-forge
recipe builds with BGEN_NANOBIND_BACKEND=private (a self-contained backend,
see python/CMakeLists.txt) and sets this same env var so the wheel metadata
doesn't falsely claim a dependency that build genuinely doesn't need.
"""

import os


def dynamic_metadata(settings, project):
    if os.environ.get("BGEN_NANOBIND_BACKEND", "shared") == "private":
        return {"dependencies": []}
    return {"dependencies": ["nanobind-backend>=1.0"]}
