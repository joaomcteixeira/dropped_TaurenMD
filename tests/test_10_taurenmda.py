# ATTENTION: tests other MATTERS because they operate over trajectory.
import os
import itertools as it
import pytest

from tauren import tauren
from tauren import load

file_path = "tests"

trajectory = os.path.join(file_path, "reference", "traj_test_PCNA.dcd")
topology = os.path.join(file_path, "reference", "topology_test.pdb")
noHOHtop = os.path.join(file_path, "reference", "mda_noHOH000.pdb")
chainid1 = os.path.join(file_path, "reference", "mda_segid_A_000.pdb")
aligned0 = os.path.join(file_path, "reference", "mda_aligned_000.pdb")
aligned50 = os.path.join(file_path, "reference", "mda_aligned_050.pdb")

def atom_predicate(line):
    
    if line.startswith("ATOM"):
        return False
    else:
        return True

def test_static_gen_pdb_name_format_1():
    s = tauren.TaurenMDAnalysis._gen_pdb_name_format(500, 'pdb')
    assert s == "{:0>3}.pdb"


def test_static_gen_pdb_name_format_2():
    s = tauren.TaurenMDAnalysis._gen_pdb_name_format(5, 'pdb')
    assert s == "{:0>1}.pdb"


def test_static_gen_selector_1():
    
    s = tauren.TaurenMDAnalysis._gen_selector(
        identifiers=[1,2,3,4],
        boolean="and",
        selection="resid",
        )
    
    assert s == "resid 1 and resid 2 and resid 3 and resid 4"

def test_static_gen_selector_2():
    
    with pytest.raises(TypeError):
        s = tauren.TaurenMDAnalysis._gen_selector(
            identifiers=5,
            boolean="and",
            selection="resid",
            )


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


def test_frameslice_1():
    
    traj.frame_slice(
        start=1,
        end=10,
        step=1,
        )
    
    assert traj.n_frames == 10
    return
    
    
def test_frames_to_slice_1():
    assert traj._gen_frame_slicer_from_string("1") == slice(0, 1, 1)


def test_frames_to_slice_2():
    assert traj._gen_frame_slicer_from_string("1:") == slice(0, 100, 1)


def test_frames_to_slice_3():
    assert traj._gen_frame_slicer_from_string(":50") == slice(0, 50, 1)


def test_frames_to_slice_4():
    assert traj._gen_frame_slicer_from_string("1:10") == slice(0, 10, 1)


def test_frames_to_slice_5():
    assert traj._gen_frame_slicer_from_string("1:50:2") == slice(0, 50, 2)


def test_frames_to_slice_6():
    assert traj._gen_frame_slicer_from_string("::2") == slice(0, 100, 2)


def test_frames_to_slice_7():
    with pytest.raises(ValueError):
        traj._gen_frame_slicer_from_string(":::")

def test_frames_to_slice_8():
    with pytest.raises(ValueError):
        traj._gen_frame_slicer_from_string("5:::")


def test_frames_to_slice_9():
    with pytest.raises(ValueError):
        traj._gen_frame_slicer_from_string(":::7")


def test_frames_to_slice_10():
    with pytest.raises(ValueError):
        traj._gen_frame_slicer_from_string("::4:7")

def test_frames_to_slice_11():
    with pytest.raises(ValueError):
        traj._gen_frame_slicer_from_string(":3:4:7")


def test_frames_to_slice_12():
    with pytest.raises(ValueError):
        traj._gen_frame_slicer_from_string("as10")


def test_frames_to_list_1():
    assert traj._get_frame_list_from_string("all") == [1,2,3,4,5,6,7,8,9,10]


def test_frames_to_list_2():
    assert traj._get_frame_list_from_string("7") == [7]


def test_frames_to_list_3():
    assert traj._get_frame_list_from_string("1,2") == [1,2]


def test_frames_to_list_4():
    assert traj._get_frame_list_from_string("1:8") == [1,2,3,4,5,6,7,8]


def test_frames_to_list_5():
    assert traj._get_frame_list_from_string("95:") == [95,96,97,98,99, 100]


def test_frames_to_list_6():
    assert traj._get_frame_list_from_string(":3") == [1,2,3]


def test_frames_to_list_6():
    assert traj._get_frame_list_from_string("1:10:2") == [1,3,5,7,9]


def test_frameslice_4():
    
    traj.frame_slice(
        start=None,
        end=None,
        step=1,
        )
    
    assert traj.n_frames == 100
    return


def test_rmv_solvent_selector_default():
    assert traj._rmv_solvent_selector == "all"


def test_default_atom_selection():
    assert traj.atom_selection == "all and all"


def test_frame_with_solvent():
    
    traj.frames2file(prefix="mda_with_solvent", frames="0")
    
    with open(topology, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    with open("mda_with_solvent000.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    os.remove("mda_with_solvent000.pdb")


def test_remove_solvent_selector():
    
    traj.remove_solvent()
    assert traj._rmv_solvent_selector == "(protein or nucleic)"


def test_remove_solvent_atom_selection():
    assert traj.atom_selection == "(protein or nucleic) and all"


def test_remove_solvent():
    
    traj.frames2file(prefix="_", frames="0")
    
    with open(noHOHtop, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    with open("_000.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    os.remove("_000.pdb")


def test_set_atom_selection():
    
    traj.set_atom_selection(selector="segid A")
    
    traj.frames2file(
        prefix="segid_A_",
        frames="0",
        )
    
    with open(chainid1, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    with open("segid_A_000.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    os.remove("segid_A_000.pdb")


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
    
    with open(noHOHtop, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    with open("all_000.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    os.remove("all_000.pdb")


def test_aligned_1():
    
    traj.align_traj(inplace=True)
    
    d = {
        "frames": "0",
        "prefix": "test_mda_aligned_",
        "ext": "pdb"
        }
    
    traj.frames2file(**d)
    
    with open(aligned0, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())

    with open("test_mda_aligned_000.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    os.remove("test_mda_aligned_000.pdb")

def test_aligned_2():
    
    d = {
        "frames": "50",
        "prefix": "test_mda_aligned_",
        "ext": "pdb"
        }
    
    traj.frames2file(**d)
    
    with open(aligned50, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    with open("test_mda_aligned_050.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    os.remove("test_mda_aligned_050.pdb")


def test_reset_rmv_solvent():
    traj.reset_rmv_solvent()
    assert traj.atom_selection == "all and all"

def test_rmsds_combined_1():
    
    key = traj.calc_rmsds_combined_chains()
    
    assert isinstance(key.datatype, str)
    assert isinstance(key.identifier, str)
    assert isinstance(key.filenaming, str)    
    
    traj.export_data(key)

def test_rmsds_combined_2():
    
    key = traj.calc_rmsds_separated_chains()
    
    assert isinstance(key.datatype, str)
    assert isinstance(key.identifier, str)
    assert isinstance(key.filenaming, str)  
    
    traj.export_data(key)
