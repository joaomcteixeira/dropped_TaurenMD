import numpy as np

def atom_predicate(line):
    
    if line.startswith("ATOM"):
        return False
    else:
        return True


def compare_csv(f1, f2):
    """
    f1 and f2 are file paths
    """
    
    with open(f1, 'r') as fh:
        l1 = fh.readlines()
    
    with open(f2, 'r') as fh:
        l2 = fh.readlines()
    
    # confirms that files have the same structure
    general = all(x[0] == y[0] for x, y in zip(l1, l2))
    
    # confirms comments are the same
    comments = all(x == y for x, y in zip(l1, l2) if x.startswith("#"))
    
    # confirms values are the same within a tolerance
    # rmsds may be different across computers.
    n1 = np.loadtxt(f1, delimiter=",", comments="#")
    n2 = np.loadtxt(f2, delimiter=",", comments="#")
    
    numbers = np.allclose(n1, n2, rtol=0.00002)
    
    assert general
    assert comments
    assert numbers
    
    return general and comments and numbers
