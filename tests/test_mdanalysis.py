import os
import filecmp
import itertools as it

from tauren import tauren
from tauren import load

file_path = "tests"

trajectory = os.path.join(file_path, "reference", "traj_test_PCNA.dcd")
topology = os.path.join(file_path, "reference", "topology_test.pdb")
noHOHtop = os.path.join(file_path, "reference", "mda_noHOH000.pdb")
chainid1 = os.path.join(file_path, "reference", "mda_segid_A_000.pdb")

def atom_predicate(line):
    
    if line.startswith("ATOM"):
        return True
    else:
        return False

def test_load_traj():

    traj = load.load_traj(
        trajectory,
        topology,
        traj_type="mdanalysis",
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

def test_rmv_solvent_selector_default():
    
    assert traj._rmv_solvent_selector == "all"

def test_default_atom_selection():
    
    assert traj.atom_selection == "all and all"

def test_frame_with_solvent():
    
    traj.frames2file(prefix="mda_with_solvent", frames="0")
    
    fh1 = open(topology, 'r')
    f1 = it.filterfalse(atom_predicate, fh1)
    
    fh2 = open("mda_with_solvent000.pdb", 'r')
    f2 = it.filterfalse(atom_predicate, fh2)
    
    assert all(x == y for x, y in zip(f1, f2))
    fh1.close()
    fh2.close()
    


def test_remove_solvent_selector():
    
    traj.remove_solvent()
    assert traj._rmv_solvent_selector == "(protein or nucleic)"

def test_remove_solvent_atom_selection():
    assert traj.atom_selection == "(protein or nucleic) and all"

def test_remove_solvent():
    
    traj.frames2file(prefix="_", frames="0")
    
    # with does not work because
    # >   assert all(x == y for x, y in zip(f1, f2))
    # E   ValueError: I/O operation on closed file.
    
    fh1 = open(noHOHtop, 'r')
    f1 = it.filterfalse(atom_predicate, fh1)
    
    fh2 = open("_000.pdb", 'r')
    f2 = it.filterfalse(atom_predicate, fh2)
    
    os.remove("_000.pdb")
    
    assert all(x == y for x, y in zip(f1, f2))
    fh1.close()
    fh2.close()

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


def test_set_atom_selection():
    
    traj.set_atom_selection(selector="segid A")
    
    traj.frames2file(
        prefix="segid_A_",
        frames="0",
        )
    
    fh1 = open(chainid1, 'r')
    f1 = it.filterfalse(atom_predicate, fh1)
    
    fh2 = open("segid_A_000.pdb", 'r')
    f2 = it.filterfalse(atom_predicate, fh2)
    
    os.remove("segid_A_000.pdb")
    
    assert all(x == y for x, y in zip(f1, f2))
    fh1.close()
    fh2.close()

def test_set_atom_selection_reset():
    """
    tests that atom_selection can be dynamically set.
    """
    
    traj.set_atom_selection(selector=None)
    
    assert traj.atom_selection == "(protein or nucleic) and all"

def test_set_atom_selection_2():

    traj.frames2file(
        prefix="all_",
        frames="0",
        )
    
    fh1 = open(noHOHtop, 'r')
    f1 = it.filterfalse(atom_predicate, fh1)
    
    fh2 = open("all_000.pdb", 'r')
    f2 = it.filterfalse(atom_predicate, fh2)
    
    os.remove("all_000.pdb")
    
    assert all(x == y for x, y in zip(f1, f2))
    fh1.close()
    fh2.close()

def test_reset_rmv_solvent():
    
    traj.reset_rmv_solvent()
    
    assert traj.atom_selection == "all and all"
