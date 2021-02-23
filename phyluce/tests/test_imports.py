"""Test various imports from packages"""

import phyluce


def test_phyluce_version():
    """Ensure we can successfully import"""
    assert phyluce.__version__
