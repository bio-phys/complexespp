import pytest
import yaml
import MDAnalysis as mda

from pycomplexes import addlinker
from pycomplexes.testing import datafiles


def test_get_linker_connections():
    top = yaml.load(
        """
    connections:
      0:
        domain-a: [r_beginning, B 1]
        domain-b: [deca_ala, B 2]
        type: flat
        params: {}
      1:
        domain-a: [deca_ala, B 9]
        domain-b: [r_end, B 10]
        type: flat
        params: {}
      2:
        domain-a: [r_beginning, B 1]
        domain-b: [r_end, B 10]
        type: gaussian
        params: {N: 8, bond-length: 3.81}
      3:
        domain-a: [r_beginning_2, B 1]
        domain-b: [deca_ala_2, B 2]
        type: flat
        params: {}
      4:
        domain-a: [deca_ala_2, B 9]
        domain-b: [r_end_2, B 10]
        type: flat
        params: {}
      5:
        domain-a: [r_beginning_2, B 1]
        domain-b: [r_end_2, B 10]
        type: gaussian
        params: {N: 8, bond-length: 3.81}
    """
    )
    lc = addlinker.get_linker_connections(top)

    assert len(lc) == 2

    l = lc[0]
    assert l.start.name == "r_beginning"
    assert l.start.bead == "B 1"
    assert l.stop.name == "r_end"
    assert l.stop.bead == "B 10"

    assert l.linker[0] == "deca_ala"
    assert l.linker[1] == "B 2"
    assert l.linker[2] == "B 9"

    l = lc[1]
    assert l.start.name == "r_beginning_2"
    assert l.start.bead == "B 1"
    assert l.stop.name == "r_end_2"
    assert l.stop.bead == "B 10"

    assert l.linker[0] == "deca_ala_2"
    assert l.linker[1] == "B 2"
    assert l.linker[2] == "B 9"


@pytest.fixture
def cplx(datafiles):
    with open(datafiles["chmp3.cplx"]) as fh:
        return yaml.load(fh)


@pytest.fixture
def sim(datafiles):
    return mda.Universe(datafiles["chmp3.pdb"])


def test_get_type_ags(cplx, sim):
    ags = addlinker.get_type_ags(cplx, sim, "rigid")
    assert len(ags) == 3
    assert ags[0].n_atoms == 132
    assert ags[1].n_atoms == 16
    assert ags[2].n_atoms == 16


def test_get_linkers(cplx, sim):
    linkers = addlinker.get_linkers(cplx, sim)

    assert len(linkers) == 2

    link = linkers[0]

    start = link.start
    assert start.segid == "A"
    assert start.resid == 140

    stop = link.stop
    assert stop.segid == "A"
    assert stop.resid == 157

    assert link.name == "link1"
    assert link.ag.n_atoms == 16
    assert link.ag.resids.tolist() == list(range(141, 141 + 16))
