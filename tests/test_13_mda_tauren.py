# tests TaurenTraj interface for MDAnalysis core
# ATTENTION: tests other MATTERS because they operate over trajectory.
import os
import itertools as it
import pytest

from tauren import tauren
from tauren import load

file_path = "tests"
rf = "reference"
mdtype = "mda"
trajtype = "mdanalysis"

trajectory = os.path.join(file_path, rf, "traj_test_PCNA.dcd")
topology = os.path.join(file_path, rf, "topology_test.pdb")

solvent0 = os.path.join(file_path, rf, "solvent_all_frame0.pdb")
solvent49 = os.path.join(file_path, rf, "solvent_all_frame49.pdb")

noHOH0 = os.path.join(file_path, rf, "noHOH_all_frame0.pdb")
noHOH49 = os.path.join(file_path, rf, "noHOH_all_frame49.pdb")

noHOH_chainA0 = os.path.join(file_path, rf, "noHOH_chainA_frame0.pdb")
noHOH_chainA49 = os.path.join(file_path, rf, "noHOH_chainA_frame49.pdb")

noHOH_aligned0 = os.path.join(file_path, rf, "noHOH_all_aligned_frame0.pdb")
noHOH_aligned49 = os.path.join(file_path, rf, "noHOH_all_aligned_frame49.pdb")


def atom_predicate(line):
    
    if line.startswith("ATOM"):
        return False
    else:
        return True


traj = load.load_traj(
        trajectory,
        topology,
        traj_type="mdanalysis",
        )


def test_report():
    
    report = """* Trajectory details:

    num of frames: 100
    num of residues: 1114
    num of atoms: 12160

    total time: 99.000

    time_step: 1.000 (usually ps)
        or 0.001 (usually ns)
*
"""
    
    assert report == traj.report()
    return


def test_frame_slice_1():
    
    traj.frame_slice(
        start=1,
        end=None,
        step=1,
        )
    
    assert traj.n_frames == 100


def test_frame_slice_2():
    
    traj.frame_slice(
        start=50,
        end=None,
        step=1,
        )
    
    assert traj.n_frames == 51


def test_frame_slice_3():
    
    traj.frame_slice(
        start=1,
        end=10,
        step=1,
        )
    
    assert traj.n_frames == 10


def test_frameslice_4():
    
    traj.frame_slice(
        start=None,
        end=None,
        step=1,
        )
    
    assert traj.n_frames == 100


def test_rmv_solvent_selector_default():
    assert traj._rmv_solvent_selector == "all"


def test_default_atom_selection():
    assert traj.atom_selection == "all and all"


def test_frame2file_with_solvent_1():
    
    traj.frames2file(prefix="with_solvent", frames="1")
    
    with open(solvent0, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    with open("with_solvent000.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    os.remove("with_solvent000.pdb")


def test_frame2file_with_solvent_50():
    
    traj.frames2file(prefix="with_solvent", frames="50")
    
    with open(solvent49, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    with open("with_solvent049.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    os.remove("with_solvent049.pdb")


def test_remove_solvent_selector():
    traj.remove_solvent()
    assert traj._rmv_solvent_selector == "(protein or nucleic)"


def test_remove_solvent_atom_selection():
    assert traj.atom_selection == "(protein or nucleic) and all"


def test_frame2file_noHOH_1():
    
    traj.frames2file(prefix="_", frames="1")
    
    with open(noHOH0, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    with open("_000.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    os.remove("_000.pdb")


def test_frame2file_noHOH_50():
    
    traj.frames2file(prefix="_", frames="50")
    
    with open(noHOH49, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    with open("_049.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    os.remove("_049.pdb")


def test_frame2file_chainA_1():
    
    traj.set_atom_selection(selector="segid A")
    
    traj.frames2file(
        prefix="chain_A_",
        frames="1",
        )
    
    with open(noHOH_chainA0, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    with open("chain_A_000.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    os.remove("chain_A_000.pdb")


def test_frame2file_chainA_50():
    
    traj.set_atom_selection(selector="segid A")
    
    traj.frames2file(
        prefix="chain_A_",
        frames="50",
        )
    
    with open(noHOH_chainA49, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    with open("chain_A_049.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    os.remove("chain_A_049.pdb")


def test_set_atom_selection_reset():
    """
    tests that atom_selection can be dynamically set.
    """
    
    traj.set_atom_selection(selector=None)
    
    assert traj.atom_selection == "(protein or nucleic) and all"


def test_frame2file_aligned_1():
    
    traj.align_traj(inplace=True)
    
    d = {
        "frames": "1",
        "prefix": "test_aligned_",
        "ext": "pdb"
        }
    
    traj.frames2file(**d)
    
    with open(noHOH_aligned0, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())

    with open("test_aligned_000.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    os.remove("test_aligned_000.pdb")


def test_frame2file_aligned_50():
    
    traj.align_traj(inplace=True)
    
    d = {
        "frames": "50",
        "prefix": "test_aligned_",
        "ext": "pdb"
        }
    
    traj.frames2file(**d)
    
    with open(noHOH_aligned49, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())

    with open("test_aligned_049.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    os.remove("test_aligned_049.pdb")


def test_undo_rmv_solvent_1():
    traj.undo_rmv_solvent()
    assert traj.atom_selection == "all and all"

def test_rmsds_combined_1():
    
    key = traj.calc_rmsds_combined_chains()
    
    assert len(traj.observables) == 1

def test_export_data_1():
    
    traj.export_data(0)
    
    
    
    # assert isinstance(key.datatype, str)
    # assert isinstance(key.identifier, str)
    # assert isinstance(key.filenaming, str)    
    
    # traj.export_data(key)

# def test_rmsds_combined_2():
    
    # key = traj.calc_rmsds_separated_chains()
    
    # assert isinstance(key.datatype, str)
    # assert isinstance(key.identifier, str)
    # assert isinstance(key.filenaming, str)  
    
    # traj.export_data(key)
