# tests MDtraj dunders
import os
import string
import pytest

from tauren import (tauren, load)

file_path = "tests"
mdtype = "mdt"
trajtype = "mdtraj"

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
    
    assert chainlist == list(range(traj.trajectory.n_chains))


def test_rmv_solvent_1():
    
    trajref = traj.trajectory
    
    assert trajref == traj.trajectory


def test_rmv_solvent_2():
    
    trajref = traj.trajectory
    
    traj.remove_solvent()
    
    assert trajref != traj.trajectory
    
    traj.undo_rmv_solvent()
    
    assert trajref == traj.trajectory


def test_frames_to_list_1():
    assert traj._get_frame_list("all") == list(range(1, 101, 1))


def test_frames2file_1():
    
    traj._frames2file([1], "__{}.pdb")
    os.remove("__0.pdb")
