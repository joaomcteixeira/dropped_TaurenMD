import pytest

from tauren.tauren import TrajObservables

def test_1():
    
    a = TrajObservables()
    a.append(name="testing", data=1)
    
    assert a == [{"name": "testing", "data": 1}]


def test_2():
    
    a = TrajObservables()
    a.append(name="testing1", data=1)
    a.append(name="testing2", data=2)
    
    ref = [{"name": "testing1", "data": 1}, {"name": "testing2", "data": 2}]
    
    assert a == ref

def test_3():
    
    a = TrajObservables()
    a.append("just a string")
    
    assert a == [{"data_1": "just a string"}]

def test_4():
    
    a = TrajObservables()
    a.append(1)
    
    assert a == [{"data_1": 1}]

def test_5():
    
    a = TrajObservables()
    a.append("a string", 1, name="testing1", data=1)
    
    assert a == [{"data_1": "a string", "data_2": 1, "name": "testing1", "data": 1}]

def test_6():
    
    a = TrajObservables()
    a.append("a string", 1, {"data_1": "just a string"})
    
    assert a == [{"data_1": "a string", "data_2": 1, "data_3": {"data_1": "just a string"}}]

def test_7():
    
    a = TrajObservables()
    with pytest.raises(ValueError):
        
        a.append("a string", data_1=1)
