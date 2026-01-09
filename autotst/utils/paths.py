from __future__ import annotations

from pathlib import Path


def project_root() -> Path:
    # Adjust if this file moves deeper/shallower in the tree.
    return Path(__file__).resolve().parents[2]


def test_dir() -> Path:
    return project_root() / "test"


def test_bin_dir() -> Path:
    return test_dir() / "bin"


def test_data_dir() -> Path:
    return test_bin_dir() / "log-files"


def database_dir() -> Path:
    return project_root() / "database"
