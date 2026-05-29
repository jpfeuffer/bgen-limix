"""S3/UPath credential helpers for cbgen."""

from __future__ import annotations

import contextlib
import os


@contextlib.contextmanager
def s3_env_from_storage_options(filepath):
    """Temporarily set AWS env vars from a UPath/fsspec ``storage_options`` dict.

    UPath objects (from ``universal_pathlib``) carry credentials as
    ``storage_options``, using the same keys as s3fs / fsspec:

    * ``key``          → ``AWS_ACCESS_KEY_ID``
    * ``secret``       → ``AWS_SECRET_ACCESS_KEY``
    * ``token``        → ``AWS_SESSION_TOKEN``
    * ``endpoint_url`` → ``AWS_ENDPOINT_URL``
    * ``region_name``  → ``AWS_DEFAULT_REGION``
    * ``anon=True``    → ``AWS_NO_SIGN_REQUEST=1``

    ``client_kwargs`` (a nested dict used by s3fs) is also checked for
    ``endpoint_url`` and ``region_name``.

    The env vars are restored to their original values (or removed) when
    the context exits, so this is safe to use concurrently as long as
    each thread/process opens its file handle within its own context.
    """
    opts = getattr(filepath, "storage_options", None)
    if not opts:
        yield
        return

    to_set: dict[str, str] = {}

    for fsspec_key, env_var in [
        ("key",          "AWS_ACCESS_KEY_ID"),
        ("secret",       "AWS_SECRET_ACCESS_KEY"),
        ("token",        "AWS_SESSION_TOKEN"),
        ("endpoint_url", "AWS_ENDPOINT_URL"),
        ("region_name",  "AWS_DEFAULT_REGION"),
    ]:
        val = opts.get(fsspec_key)
        if val:
            to_set[env_var] = str(val)

    # s3fs also accepts endpoint_url / region_name nested inside client_kwargs
    for fsspec_key, env_var in [
        ("endpoint_url", "AWS_ENDPOINT_URL"),
        ("region_name",  "AWS_DEFAULT_REGION"),
    ]:
        val = opts.get("client_kwargs", {}).get(fsspec_key)
        if val and env_var not in to_set:
            to_set[env_var] = str(val)

    if opts.get("anon"):
        to_set["AWS_NO_SIGN_REQUEST"] = "1"

    saved = {k: os.environ.get(k) for k in to_set}
    os.environ.update(to_set)
    try:
        yield
    finally:
        for k, old in saved.items():
            if old is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = old
