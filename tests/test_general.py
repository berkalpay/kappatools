import pytest

from kappatools import kasnap


@pytest.fixture
def example_mixture():
    return kasnap.SnapShot("tests/data/large.ka")


@pytest.fixture
def example_complex(example_mixture):
    return example_mixture.complexes[0]


def test_complex_size(example_complex):
    assert example_complex.size == 43