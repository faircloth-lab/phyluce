"""Test various imports from packages"""


def test_phyluce_version():
    """Ensure we can successfully import"""
    import phyluce
    assert phyluce.__version__.startswith("git")
