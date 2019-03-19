# ATTENTION: tests other MATTERS because they operate over trajectory.
import os
import filecmp
import itertools as it

from tauren import tauren
from tauren import load

file_path = "tests"

trajectory = os.path.join(file_path, "reference", "traj_test_PCNA.dcd")
topology = os.path.join(file_path, "reference", "topology_test.pdb")
noHOHtop = os.path.join(file_path, "reference", "mdt_noHOH000.pdb")
chainid1 = os.path.join(file_path, "reference", "mdt_chainid_0_000.pdb")

def atom_predicate(line):
    
    if line.startswith("ATOM"):
        return False
    else:
        return True

def test_load_traj():

    traj = load.load_traj(
        trajectory,
        topology,
        traj_type="mdtraj",
        )

    return traj

traj = test_load_traj()

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

def test_remove_solvent():
    
    traj.remove_solvent()
    
    traj.frames2file(prefix="_", frames="0")
    
    with open(noHOHtop, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    with open("_000.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    os.remove("_000.pdb")


def test_frameslice_1():
    
    traj.frame_slice(
        start=1,
        end=10,
        step=1,
        )
    
    assert traj.n_frames == 10
    return

def test_frameslice_2():
    
    traj.frame_slice(
        start=1,
        end=None,
        step=1,
        )
    
    assert traj.n_frames == 100
    return

def test_frameslice_3():
    
    traj.frame_slice(
        start=50,
        end=None,
        step=1,
        )
    
    assert traj.n_frames == 51
    return

def test_frameslice_4():
    
    traj.frame_slice(
        start=None,
        end=None,
        step=1,
        )
    
    assert traj.n_frames == 100
    return

def test_default_atom_selection():
    
    assert traj.atom_selection == "all"


def test_set_atom_selection():
    
    traj.set_atom_selection(selector="chainid 0")
    
    traj.frames2file(
        prefix="chainid_0_",
        frames="0",
        )
    
    with open(chainid1, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    with open("chainid_0_000.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    os.remove("chainid_0_000.pdb")


def test_set_atom_selection_2():
    """
    tests that atom_selection can be dynamically set.
    """
    
    traj.set_atom_selection(selector=None)
    
    traj.frames2file(
        prefix="all_",
        frames="0",
        )
    
    with open(noHOHtop, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    with open("all_000.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    os.remove("all_000.pdb")

def test_rmsds_combined_1():
    
    traj.calc_rmsds_combined_chains()
    
    
