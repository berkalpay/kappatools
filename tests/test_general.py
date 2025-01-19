import pytest

from kappatools import kasnap, kamol, kamatch, kagraph


@pytest.fixture
def example_mixture():
    return kasnap.SnapShot("tests/data/large.ka")


@pytest.fixture
def example_complex(example_mixture):
    return example_mixture.complexes[0]


def test_complex_size(example_complex):
    assert example_complex.size == 43


def test_decode_isomorphic(example_complex):
    decoded_data = kamol.Kappa().decode(
        example_complex.canonical, example_complex.system_views
    )
    decoded_complex = kamol.KappaComplex(decoded_data)
    assert kamatch.SiteGraphMatcher().isomorphic(example_complex, decoded_complex)


def test_cycle():
    complex = kamol.KappaComplex(
        "A(r[6] l[1]),A(r[1] l[2]),A(r[2] l[3]),A(r[3] l[4]),A(r[4] l[5]),A(r[5] l[6])"
    )
    cycle = kagraph.KappaGraph(complex).get_cycle()
    assert len(cycle) == 6


def test_kamol():
    # parser
    kappa = kamol.Kappa()
    # a simple Kappa string
    s1 = " A(o[1], p[2] t{p}[3]), B(x[1] y[2] z[.]), C(w[3], y[z.B])"
    agents = kappa.parser(s1)
    c = kamol.KappaExpression(agents)
    out = c.kappa_expression()
    print(f"expression:\n{out}")
    out = c.show(internal=True)
    print(f"internal representation:\n{out}")

    print("--------------")

    s2 = " x222667:P(a1[.] a2[.] a3[.] d[1]), x251065:P(a1[.] a2[.] a3[.] d[1])"
    agents = kappa.parser(s2)
    c = kamol.KappaMolecule(agents)
    out = c.kappa_expression()
    print(f"expression:\n{out}")
    out = c.show(internal=True)
    print(f"internal representation:\n{out}")

    print("--------------")

    print("Reading snapshot")
    snapshot = kasnap.SnapShot(file="tests/data/1784.ka")
    print("Done reading")
    print(snapshot.snap_report())

    # print("--------------")
    #
    SGM = kamatch.SiteGraphMatcher()
    #
    # for molecule in snapshot.complexes:
    #     canonical = molecule.canonical
    #     print(canonical)
    #     molecule_ = KappaComplex(kappa.decode(canonical, snapshot.local_views))
    #     print(f'decoded is isomorphic to original: {SGM.isomorphic(molecule, molecule_)}')

    print("--------------")

    print("Reading snapshot")
    snapshot = kasnap.SnapShot(file="tests/data/1773.ka")
    print("Done reading")
    print(snapshot.snap_report())

    example = snapshot.complexes[0]
    print(example.show(label=True, wrap=100))
    print("same thing")
    s = kamol.KappaMolecule(example.agents)
    print(s.show(label=True, wrap=100))

    print("> remapped______________")
    example.remap(change="randomize", id_shift=100)
    print(example.show(label=True, wrap=100))

    s = kamol.KappaMolecule(example.agents)
    print("agent ids have been normalized")
    print(s.show(label=True, wrap=100))

    print(s.summary(internal=True, show_bonds=True, reactivity=True))

    # import kasig as sig

    data = "A(a1[1] a2[2] a3[3] c[8]), A(a1[1] a2[2] a3[3] c[4]) A(a1[5] a2[6] a3[7] c[4]), A(a1[5] a2[6] a3[7] c[8])"
    x1 = kamol.KappaComplex(data)
    print(x1.kappa_expression())
    print(x1.canonical)
    out = kappa.decode(x1.canonical, x1.local_view_index)
    print(out)
    print(
        f"decoded is isomorphic to original: {SGM.isomorphic(kamol.KappaComplex(out), x1)}"
    )
    expression = kamol.Canonical2Expression(x1.canonical, x1.local_view_index)
    print(expression)
    print(f"decoded is isomorphic to original: {SGM.isomorphic(expression, x1)}")
