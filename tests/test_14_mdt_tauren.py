# tests TaurenTraj interface for MDAnalysis core
# ATTENTION: tests other MATTERS because they operate over trajectory.
import os
import itertools as it
import pytest
import numpy as np
import json
import copy

from tauren import tauren
from tauren import load

file_path = "tests"
rf = "reference"
mdtype = "mdtraj"
trajtype = "mdtraj"

trajectory = os.path.join(file_path, rf, "traj_test_PCNA.dcd")
topology = os.path.join(file_path, rf, "topology_test.pdb")

solvent0 = os.path.join(file_path, rf, mdtype, "solvent_all_frame0.pdb")
solvent49 = os.path.join(file_path, rf, mdtype, "solvent_all_frame49.pdb")

noHOH0 = os.path.join(file_path, rf, mdtype, "noHOH_all_frame0.pdb")
noHOH49 = os.path.join(file_path, rf, mdtype, "noHOH_all_frame49.pdb")

noHOH_chainA0 = os.path.join(file_path, rf, mdtype, "noHOH_chainA_frame0.pdb")
noHOH_chainA49 = os.path.join(
    file_path,
    rf,
    mdtype,
    "noHOH_chainA_frame49.pdb",
    )

# noHOH_aligned0 = os.path.join(
    # file_path,
    # rf,
    # mdtype,
    # "noHOH_all_aligned_frame0.pdb",
    # )
# noHOH_aligned49 = os.path.join(
    # file_path,
    # rf,
    # mdtype,
    # "noHOH_all_aligned_frame49.pdb",
    # )

rmsd_solvent = os.path.join(
    file_path,
    rf,
    mdtype,
    "rmsds_combined_chains_-all-_and_-all-_all.csv"
    )

rmsd_noHOH = os.path.join(
    file_path,
    rf,
    mdtype,
    "rmsds_combined_chains_-protein_or_nucleic-_and_-all-_all.csv"
    )

rmsd_noHOH_A = os.path.join(
    file_path,
    rf,
    mdtype,
    "rmsds_combined_chains_-protein_or_nucleic-_and_-all-_A.csv"
    )

rmsd_sep_solvent = os.path.join(
    file_path,
    rf,
    mdtype,
    "rmsds_separated_chains_-all-_and_-all-_all.csv"
    )

rmsd_sep_noHOH = os.path.join(
    file_path,
    rf,
    mdtype,
    "rmsds_separated_chains_-protein_or_nucleic-_and_-all-_all.csv"
    )

rmsd_sep_noHOH_A = os.path.join(
    file_path,
    rf,
    mdtype,
    "rmsds_separated_chains_-protein_or_nucleic-_and_-all-_A.csv"
    )

dataindex6 = os.path.join(file_path, rf, "data_index_6.json")


def atom_predicate(line):
    
    if line.startswith("ATOM"):
        return False
    else:
        return True


traj = load.load_traj(
        trajectory,
        topology,
        traj_type=trajtype,
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
    with pytest.raises(AttributeError):
        traj._rmv_solvent_selector


def test_default_atom_selection():
    assert traj.atom_selection == "(all)"


def test_frame2file_with_solvent_1():
    
    traj.frames2file(prefix="with_solvent", frames="1")
    
    with open(solvent0, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    with open("with_solvent000.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    # os.remove("with_solvent000.pdb")


def test_frame2file_with_solvent_50():
    
    traj.frames2file(prefix="with_solvent", frames="50")
    
    with open(solvent49, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    with open("with_solvent049.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    # os.remove("with_solvent049.pdb")


def test_remove_solvent_selector():
    # t1 = copy.copy(traj.trajectory)
    traj.remove_solvent()
    # t2 = copy.copy(traj.trajectory)
    # assert t1 == t2


def test_remove_solvent_atom_selection():
    assert traj.atom_selection == "(all)"

def test_frame2file_noHOH_1():
    
    traj.frames2file(prefix="_", frames="1")
    
    with open(noHOH0, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    with open("_000.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    # os.remove("_000.pdb")


def test_frame2file_noHOH_50():
    
    traj.frames2file(prefix="_", frames="50")
    
    with open(noHOH49, 'r') as fh:
        f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    with open("_049.pdb", 'r') as fh:
        f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    assert all(x == y for x, y in zip(f1, f2))
    # os.remove("_049.pdb")


# def test_chainid0():
    
    # traj.set_atom_selection(selector="chainid 0")
    # assert traj.atom_selection == "(chainid 0)"


# def test_frame2file_chainA_1():
    
    # # traj.set_atom_selection(selector="chainid 0")
    
    # traj.frames2file(
        # prefix="chain_A_",
        # frames="1",
        # )
    
    # with open(noHOH_chainA0, 'r') as fh:
        # f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    # with open("chain_A_000.pdb", 'r') as fh:
        # f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    # assert all(x == y for x, y in zip(f1, f2))
    # os.remove("chain_A_000.pdb")


# def test_frame2file_chainA_50():
    
    # # traj.set_atom_selection(selector="chainid 0")
    
    # traj.frames2file(
        # prefix="chain_A_",
        # frames="50",
        # )
    
    # with open(noHOH_chainA49, 'r') as fh:
        # f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    # with open("chain_A_049.pdb", 'r') as fh:
        # f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    # assert all(x == y for x, y in zip(f1, f2))
    # os.remove("chain_A_049.pdb")


# def test_set_atom_selection_reset():
    # """
    # tests that atom_selection can be dynamically set.
    # """
    
    # traj.set_atom_selection(selector=None)
    
    # assert traj.atom_selection == "(all)"


# def test_frame2file_aligned_1():
    
    # assert "not implemented" == traj.align_traj(inplace=True)
    
    # d = {
        # "frames": "1",
        # "prefix": "test_aligned_",
        # "ext": "pdb"
        # }
    
    # traj.frames2file(**d)
    
    # with open(noHOH_aligned0, 'r') as fh:
        # f1 = it.filterfalse(atom_predicate, fh.readlines())

    # with open("test_aligned_000.pdb", 'r') as fh:
        # f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    # assert all(x == y for x, y in zip(f1, f2))
    # os.remove("test_aligned_000.pdb")


# def test_frame2file_aligned_50():
    
    # traj.align_traj(inplace=True)
    
    # d = {
        # "frames": "50",
        # "prefix": "test_aligned_",
        # "ext": "pdb"
        # }
    
    # traj.frames2file(**d)
    
    # with open(noHOH_aligned49, 'r') as fh:
        # f1 = it.filterfalse(atom_predicate, fh.readlines())

    # with open("test_aligned_049.pdb", 'r') as fh:
        # f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    # assert all(x == y for x, y in zip(f1, f2))
    # os.remove("test_aligned_049.pdb")


# def test_undo_rmv_solvent_1():
    # traj.undo_rmv_solvent()
    # # assert traj.atom_selection == "(all) and (all)"

# def test_frame2file_noHOH_1_2():
    
    # traj.frames2file(prefix="_", frames="1")
    
    # with open(noHOH0, 'r') as fh:
        # f1 = it.filterfalse(atom_predicate, fh.readlines())
    
    # with open("_000.pdb", 'r') as fh:
        # f2 = it.filterfalse(atom_predicate, fh.readlines())
    
    # assert all(x == y for x, y in zip(f1, f2))
    # os.remove("_000.pdb")


# def test_cross_data_1():
    
    # a0 = np.loadtxt(rmsd_solvent, delimiter=",")
    # a1 = np.loadtxt(rmsd_noHOH, delimiter=",")
    # a2 = np.loadtxt(rmsd_noHOH_A, delimiter=",")
    
    # assert np.array_equal(a0[:, 0], a1[:, 0])
    # assert np.array_equal(a0[:, 0], a2[:, 0])
    
    # assert not(np.array_equal(a0[:, 1], a1[:, 1]))
    # assert not(np.array_equal(a0[:, 1], a2[:, 1]))
    # assert not(np.array_equal(a1[:, 1], a2[:, 1]))


# def test_cross_data_2():
    
    # a0 = np.loadtxt(rmsd_sep_solvent, delimiter=",")
    # a1 = np.loadtxt(rmsd_sep_noHOH, delimiter=",")
    # a2 = np.loadtxt(rmsd_sep_noHOH_A, delimiter=",")
    
    # assert np.array_equal(a0[:, 0], a1[:, 0])
    # assert np.array_equal(a0[:, 0], a2[:, 0])
    
    # assert np.array_equal(a0[:, 1], a1[:, 1])
    # assert np.array_equal(a0[:, 1], a2[:, 1])
    
    # assert np.array_equal(a0[:, 2], a1[:, 2])
    # assert np.array_equal(a0[:, 3], a1[:, 3])
    
    # assert a0.shape[1] > a1.shape[1] > a2.shape[1]


# def test_cross_data_3():
    
    # a0 = np.loadtxt(rmsd_noHOH_A, delimiter=",")
    # a1 = np.loadtxt(rmsd_sep_noHOH_A, delimiter=",")
    
    # assert np.array_equal(a0[:, 1], a1[:, 1])


# def test_rmsds_combined_1():
    
    # key = traj.calc_rmsds_combined_chains()
    
    # assert len(traj.observables) == 1


# def test_export_data_1():
    
    # traj.export_data(0)
    
    # with open(rmsd_solvent, 'r') as fh:
        # f1 = fh.readlines()
    
    # with open("rmsds_combined_chains_-all-_and_-all-_all.csv", 'r') as fh:
        # f2 = fh.readlines()
    
    # assert f1 == f2
    # os.remove("rmsds_combined_chains_-all-_and_-all-_all.csv")


# def test_rmsds_combined_2():
    
    # traj.remove_solvent()
    
    # key = traj.calc_rmsds_combined_chains()
    
    # assert len(traj.observables) == 2


# def test_export_data_2():
    
    # traj.export_data(1)
    
    # with open(rmsd_noHOH, 'r') as fh:
        # f1 = fh.readlines()
    
    # rmsdtest = "rmsds_combined_chains_-protein_or_nucleic-_and_-all-_all.csv"
    
    # with open(rmsdtest, 'r') as fh:
        # f2 = fh.readlines()
    
    # assert f1 == f2
    # os.remove(rmsdtest)


# def test_rmsds_combined_3():
    
    # key = traj.calc_rmsds_combined_chains(chains="0")
    
    # assert len(traj.observables) == 3


# def test_export_data_3():
    
    # traj.export_data(2)
    
    # with open(rmsd_noHOH_A, 'r') as fh:
        # f1 = fh.readlines()
    
    # rmsdtest = "rmsds_combined_chains_-protein_or_nucleic-_and_-all-_A.csv"
    
    # with open(rmsdtest, 'r') as fh:
        # f2 = fh.readlines()
    
    # assert f1 == f2
    # os.remove(rmsdtest)


# def test_rmsds_separated_1():
    
    # traj.undo_rmv_solvent()
    
    # key = traj.calc_rmsds_separated_chains()
    
    # assert len(traj.observables) == 4


# def test_export_data_4():
    
    # traj.export_data(3)
    
    # with open(rmsd_sep_solvent, 'r') as fh:
        # f1 = fh.readlines()
    
    # with open("rmsds_separated_chains_-all-_and_-all-_all.csv", 'r') as fh:
        # f2 = fh.readlines()
    
    # assert f1 == f2
    # os.remove("rmsds_separated_chains_-all-_and_-all-_all.csv")


# def test_rmsds_separated_2():
    
    # traj.remove_solvent()
    
    # key = traj.calc_rmsds_separated_chains()
    
    # assert len(traj.observables) == 5


# def test_export_data_5():
    
    # traj.export_data(4)
    
    # with open(rmsd_sep_noHOH, 'r') as fh:
        # f1 = fh.readlines()
    
    # rmsdtest = "rmsds_separated_chains_-protein_or_nucleic-_and_-all-_all.csv"
    
    # with open(rmsdtest, 'r') as fh:
        # f2 = fh.readlines()
    
    # assert f1 == f2
    # # os.remove(rmsdtest)


# def test_rmsds_separated_3():
    
    # key = traj.calc_rmsds_separated_chains(chains="0")
    
    # assert len(traj.observables) == 6


# def test_export_data_6():
    
    # traj.export_data(5)
    
    # with open(rmsd_sep_noHOH_A, 'r') as fh:
        # f1 = fh.readlines()
    
    # rmsdtest = "rmsds_separated_chains_-protein_or_nucleic-_and_-all-_A.csv"
    
    # with open(rmsdtest, 'r') as fh:
        # f2 = fh.readlines()
    
    # assert f1 == f2
    # os.remove(rmsdtest)

# def test_gen_export_file_name_1():
    
    # index = 0
    # file_name = None
    # prefix = None
    # suffix = None
    
    # name = traj._gen_export_file_name(
        # index,
        # file_name=file_name,
        # prefix=prefix,
        # suffix=suffix,
        # )

    # assert name == "rmsds_combined_chains_-all-_and_-all-_all.csv"


# def test_gen_export_file_name_2():
    
    # index = 0
    # file_name = None
    # prefix = "yeahh"
    # suffix = None
    
    # name = traj._gen_export_file_name(
        # index,
        # file_name=file_name,
        # prefix=prefix,
        # suffix=suffix,
        # )
    
    # assert name == "yeahh_rmsds_combined_chains_-all-_and_-all-_all.csv"


# def test_gen_export_file_name_3():
    
    # index = 0
    # file_name = None
    # prefix = "yeahh"
    # suffix = "txt"
    
    # name = traj._gen_export_file_name(
        # index,
        # file_name=file_name,
        # prefix=prefix,
        # suffix=suffix,
        # )
    
    # assert name == "yeahh_rmsds_combined_chains_-all-_and_-all-_all.txt"


# def test_gen_export_file_name_4():
    
    # index = 0
    # file_name = "superfile"
    # prefix = "yeahh"
    # suffix = "txt"
    
    # name = traj._gen_export_file_name(
        # index,
        # file_name=file_name,
        # prefix=prefix,
        # suffix=suffix,
        # )
    
    # assert name == "superfile.txt"


# def test_export_json_1():
    
    # traj.observables.append(
        # "this is super data",
        # info1 = "with some info",
        # array1 = np.arange(20),
        # )
    
    # traj.export_data(
        # 6,
        # tojson=True,
        # suffix="json",
        # header="""
        # Super header:
        # I am not free for doing what I want,
        # I am free for loving  what I do,
        # - someone unknown -
        # """,
        # )
    
    # with open("data_index_6.json", 'r') as fh:
        # ind6 = json.load(fh)
    
    # with open(dataindex6, 'r') as fh:
        # refjson6 = json.load(fh)
    
    # assert ind6 == refjson6
    # os.remove("data_index_6.json")
