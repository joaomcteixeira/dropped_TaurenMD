import pytest

from tauren.tauren import TrajObservables


def test_key_gen_1():
    
    key = TrajObservables.gen_key(
        "separated_RMSDs",
        "resid 1:50",
        ["A", "B"],
        )
    
    key_ref = TrajObservables.StorageKey(
        datatype="separated_RMSDs",
        identifier="resid-1-50_A,resid-1-50_B",
        filenaming="separated_RMSDs_resid-1-50_for_A-B",
        )
        
    assert key_ref == key

def test_key_gen_2():
    
    key = TrajObservables.gen_key(
        "combined_RMSDs",
        "resid 1:50",
        "A,B,C",
        )
    
    key_ref = TrajObservables.StorageKey(
        datatype="combined_RMSDs",
        identifier="resid-1-50_A-B-C",
        filenaming="combined_RMSDs_resid-1-50_for_A-B-C",
        )
        
    assert key_ref == key

def test_key_gen_3():
    
    with pytest.raises(TypeError):
        TrajObservables.gen_key(
            ["combined_RMSDs"],
            "resid 1:50",
            "A-B-C",
            )

def test_key_gen_4():
    
    with pytest.raises(TypeError):
        TrajObservables.gen_key(
            "combined_RMSDs",
            ["resid 1:50"],
            "A-B-C",
            )

def test_key_gen_5():
    
    with pytest.raises(TypeError):
        TrajObservables.gen_key(
            "combined_RMSDs",
            "resid 1:50",
            1,
            )

def test_store_data_1():
    
    key = TrajObservables.StorageKey(
        datatype="dummy_data",
        identifier="dddumy",
        filenaming="plot_dummy",
        )
    
    data = TrajObservables.StorageData(
        data="I am just a string",
        columns="nothing here",
        )
    
    to = TrajObservables()
    
    to.store(key, data)
    
def test_store_data_2():
    
    key = TrajObservables.StorageKey(
        datatype="dummy_data",
        identifier="dddumy",
        filenaming="plot_dummy",
        )
    
    data = "I am just a string"
    
    to = TrajObservables()
    
    with pytest.raises(TypeError):
        to.store(key, data)

def test_store_data_3():
    
    key = "I am just a string"
    
    data = TrajObservables.StorageData(
        data="I am just a string",
        columns="nothing here",
        )
    
    
    to = TrajObservables()
    
    with pytest.raises(TypeError):
        to.store(key, data)
