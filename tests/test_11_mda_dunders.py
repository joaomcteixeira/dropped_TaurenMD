# tests MDAnalysis dunders
import os
import string
import pytest

from tauren import (tauren, load)

file_path = "tests"
mdtype = "mda"
trajtype = "mdanalysis"

trajectory = os.path.join(file_path, "reference", "traj_test_PCNA.dcd")
topology = os.path.join(file_path, "reference", "topology_test.pdb")
# noHOHtop = os.path.join(file_path, "reference", mdtype, "noHOH000.pdb")
# chainid1 = os.path.join(file_path, "reference", mdtype, "segid_A_000.pdb")
# aligned0 = os.path.join(file_path, "reference", mdtype, "aligned_000.pdb")
# aligned50 = os.path.join(file_path, "reference", mdtype, "aligned_050.pdb")

traj = load.load_traj(
    trajectory,
    topology,
    traj_type=trajtype,
    )


def test_gen_chain_list_all():
    
    chainlist = traj._gen_chain_list("all")
    
    assert chainlist == list(string.ascii_letters + string.digits)


def test_rmv_solvent_1():
    
    traj.remove_solvent()
    
    assert traj.atom_selection == "(protein or nucleic) and all"


def test_rmv_solvent_2():
    
    traj.remove_solvent(exclude=["Cl"])
    
    assert traj.atom_selection == "(protein or nucleic or name Cl) and all"


def test_undo_rmv_solvent_1():
    
    traj.undo_rmv_solvent()
    
    assert traj.atom_selection == "all and all"


def test_set_atom_selection_1():
    
    traj.set_atom_selection("resid 1")
    
    assert traj.atom_selection == "all and (resid 1)"


def test_set_atom_selection_2():
    
    traj.set_atom_selection("resid 1 or resid 50")
    
    assert traj.atom_selection == "all and (resid 1 or resid 50)"


def test_rmvsol_atom_1():
    
    traj.remove_solvent()
    
    assert traj.atom_selection == "(protein or nucleic) and (resid 1 or resid 50)"


def test_revert_atom_selection():
    
    traj.set_atom_selection(None)
    
    # assert traj._rmv_solvent_selector == "protein or nucleic"
    assert traj.atom_selection == "(protein or nucleic) and all"


def test_image_molecules_1():
    
    assert traj.image_molecules() == "not implemented"


def test_filter_existent_selectors_1():
    
    sel = ["segid A", "segid B", "segid Z"]
    
    result = traj._filter_existent_selectors(sel)
    
    assert result == ["segid A", "segid B"]

def test_filter_existent_selectors_3():
    
    sel = ["segid A", "segid B", "segid Z"]
    
    traj.set_atom_selection("segid A")
    result = traj._filter_existent_selectors(sel)
    
    assert result == ["segid A"]

def test_filter_existent_selectors_4():
    
    sel = ["segid A", "segid B", "segid Z"]
    traj.set_atom_selection("all")
    result = traj._filter_existent_selectors(sel)
    
    assert result == ["segid A", "segid B"]


def test_filter_existent_selectors_2():
    
    sel = "segid A, segid B, segid Z"
    
    with pytest.raises(TypeError):
        traj._filter_existent_selectors(sel)


def test_frames_to_list_1():
    assert traj._get_frame_list("all") == list(range(1, 101, 1))
    

def test_frames2file_1():
    
    traj._frames2file([1], "__{}.pdb")
    os.remove("__0.pdb")
