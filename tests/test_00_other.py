
# coding: utf-8

# In[1]:


import os
import itertools as it


# In[2]:


def atom_predicate(line):
    
    if line.startswith("ATOM"):
        return False
    else:
        return True


# In[3]:


file1 = os.path.join("tests", "reference", "file1.txt")
file2 = os.path.join("tests", "reference", "file2.txt")
file3 = os.path.join("tests", "reference", "file3.txt")
with open(file1, 'r') as fh:
    fh1 = fh.readlines()

with open(file2, 'r') as fh:
    fh2 = fh.readlines()

with open(file3, 'r') as fh:
    fh3 = fh.readlines()


# In[4]:


def test_files_equal_ATOM():
    
    f1 = it.filterfalse(atom_predicate, fh1)
    f2 = it.filterfalse(atom_predicate, fh2)
    
    assert all(x == y for x, y in zip(f1, f2))


# In[5]:


def test_files_not_equal_ATOM():
    
    f1 = it.filterfalse(atom_predicate, fh1)
    f3 = it.filterfalse(atom_predicate, fh3)
    
    assert not(all(x == y for x, y in zip(f1, f3)))    

