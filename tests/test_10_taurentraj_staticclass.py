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

def test_gen_chain_list_11():
    
    chainlist = traj._gen_chain_list("A,,B,C,")
    
    assert chainlist == ["A", "B", "C"]


def test_frames_to_list_1():
    assert traj._get_frame_list("7") == [7]

def test_frames_to_list_2():
    assert traj._get_frame_list(7) == [7]


def test_frames_to_list_3():
    assert traj._get_frame_list("1,2") == [1,2]

# tests both frames_to_list and _gen_frame_slices_from_string
def test_frames_to_list_4():
    assert traj._get_frame_list("1:8") == [1,2,3,4,5,6,7,8]


def test_frames_to_list_5():
    assert traj._get_frame_list("95:") == [95,96,97,98,99, 100]


def test_frames_to_list_6():
    assert traj._get_frame_list(":3") == [1,2,3]


def test_frames_to_list_6():
    assert traj._get_frame_list("1:10:2") == [1,3,5,7,9]


def test_frames_to_list_7():
    with pytest.raises(ValueError):
        traj._get_frame_list("1:10~")


def test_frames_to_list_8():
    with pytest.raises(ValueError):
        traj._get_frame_list("1,B")


def test_frames_to_list_9():
    with pytest.raises(TypeError):
        traj._get_frame_list(["1,B"])


def test_frames_to_list_10():
    with pytest.raises(TypeError):
        traj._get_frame_list(1.0)


def test_frames_to_list_11():
    with pytest.raises(TypeError):
        traj._get_frame_list({"foo": "bar"})


def test_check_correct_slice_1():
    
    with pytest.raises(TypeError):
        traj._check_correct_slice(0, 10, "1")


def test_check_correct_slice_2():
    
    with pytest.raises(TypeError):
        traj._check_correct_slice(0, "10", 1)


def test_check_correct_slice_3():
    
    with pytest.raises(TypeError):
        traj._check_correct_slice("0", 10, 1)


def test_check_correct_slice_4():
    
    with pytest.raises(ValueError):
        traj._check_correct_slice(0, 10, 0)


def test_check_correct_slice_5():
    
    with pytest.raises(ValueError):
        traj._check_correct_slice(10, 0, 1)


def test_check_correct_slice_6():
    
    with pytest.raises(ValueError):
        traj._check_correct_slice(0, 10, -1)


def test_gen_selector_1():
    
    result = traj._gen_selector(["A", "B"])
    
    assert result == "segid A or segid B"


def test_gen_selector_2():
    
    result = traj._gen_selector(["A", "B"], selection="resid", boolean="and")
    
    assert result == "resid A and resid B"
