from tauren import plot

def test_get_name_1():
    
    name = plot._get_fig_name("myfig", "pdf")
    
    assert name == "myfig.pdf"

def test_get_name_2():
    
    name = plot._get_fig_name("myfig.svg", "pdf")
    
    assert name == "myfig.svg"

def test_get_name_3():
    
    name = plot._get_fig_name("myfig.aaa", "pdf")
    
    assert name == "myfig.aaa"
