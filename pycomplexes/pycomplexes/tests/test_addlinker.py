import pytest
import yaml
import MDAnalysis as mda

from pycomplexes import addlinker
from pycomplexes.testing import datafiles


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
