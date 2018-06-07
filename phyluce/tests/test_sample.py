# content of test_import.py
import phyluce


def test_version():
    """Ensure we can successfully import"""
    assert phyluce.__version__.startswith("git")
