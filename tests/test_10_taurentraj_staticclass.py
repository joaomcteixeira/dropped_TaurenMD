from tauren import tauren, load


def test_static_gen_pdb_name_format_1():
    s = tauren.TaurenMDAnalysis._gen_pdb_name_format(500, 'pdb')
    assert s == "{:0>3}.pdb"


def test_static_gen_pdb_name_format_2():
    s = tauren.TaurenMDAnalysis._gen_pdb_name_format(5, 'pdb')
    assert s == "{:0>1}.pdb"


def test_static_check_chains_1():
    
    chains = ["A", "B", "C"]
    
    tauren.TaurenTraj._check_chains_argument(chains)
    assert result is None

def test_static_check_chains_2():
    
    chains = "A,B,C"
    
    result = tauren.TaurenTraj._check_chains_argument(chains)
    assert result is None

def test_static_check_chains_3():
    
    chains = "1,2,3"
    
    result = tauren.TaurenTraj._check_chains_argument(chains)
    assert result is None

def test_static_check_chains_4():
    
    chains = [1, 2, 3]
    
    result = tauren.TaurenTraj._check_chains_argument(chains)
    assert result is None

def test_static_check_chains_5():
    
    chains = [1, 2, "3"]
    
    result = tauren.TaurenTraj._check_chains_argument(chains)
    assert result is None

def test_static_check_chains_6():
    
    chains = 1
    
    with pytest.raises(TypeError):
        tauren.TrajObservables._check_chains_argument(chains)

def test_static_check_chains_7():
    
    chains = {"foo": "bar"}
    
    with pytest.raises(TypeError):
        tauren.TrajObservables._check_chains_argument(chains)

def test_static_check_chains_8():
    
    chains = [{"foo": "bar"}, "A"]
    
    with pytest.raises(ValueError):
        tauren.TrajObservables._check_chains_argument(chains)

def test_static_check_chains_9():
    
    chains = [[1], "A"]
    
    with pytest.raises(ValueError):
        tauren.TrajObservables._check_chains_argument(chains)
