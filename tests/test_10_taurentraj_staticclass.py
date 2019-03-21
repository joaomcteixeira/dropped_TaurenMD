import os
import pytest

from tauren import tauren, load

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


def test_static_gen_pdb_name_format_1():
    s = tauren.TaurenMDAnalysis._gen_pdb_name_format(500, 'pdb')
    assert s == "{:0>3}.pdb"


def test_static_gen_pdb_name_format_2():
    s = tauren.TaurenMDAnalysis._gen_pdb_name_format(5, 'pdb')
    assert s == "{:0>1}.pdb"


def test_gen_chain_list_1():
    
    chainlist = traj._gen_chain_list("A")
    
    assert chainlist == ["A"]


def test_gen_chain_list_2():
    
    chainlist = traj._gen_chain_list("A,B,C")
    
    assert chainlist == ["A", "B", "C"]


def test_gen_chain_list_3():
    
    chainlist = traj._gen_chain_list("A,B,1")
    
    assert chainlist == ["A", "B", "1"]


def test_gen_chain_list_4():
    
    chainlist = traj._gen_chain_list(["A", "B", "1"])
    
    assert chainlist == ["A", "B", "1"]


def test_gen_chain_list_5():
    
    chainlist = traj._gen_chain_list(["A", "B", 1])
    
    assert chainlist == ["A", "B", "1"]


def test_gen_chain_list_6():
    
    chainlist = traj._gen_chain_list(1)
    
    assert chainlist == ["1"]


def test_gen_chain_list_7():
    
    chains = {"foo": "bar"}
    
    with pytest.raises(TypeError):
        traj._gen_chain_list(chains)


def test_gen_chain_list_8():
    
    chains = ("foo", "bar")
    
    with pytest.raises(TypeError):
        traj._gen_chain_list(chains)
        

def test_gen_chain_list_9():
    
    chains = [{"foo": "bar"}, "A"]
    
    with pytest.raises(ValueError):
        traj._gen_chain_list(chains)


def test_gen_chain_list_10():
    
    chains = [[1], "A"]
    
    with pytest.raises(ValueError):
        traj._gen_chain_list(chains)
