from __future__ import annotations

import importlib.metadata

import cats as m


def test_version():
    assert importlib.metadata.version("cats") == m.__version__
