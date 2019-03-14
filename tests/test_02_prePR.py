import os

from tauren import version

with open("README.rst", 'r') as fh:
    lines = fh.readlines()
    readme = lines[8].split("&color")[0][80:].strip()

with open(os.path.join("install", "messages.py")) as fh:
    for line in fh: 
        if line.startswith("# http://patorjk.com"):
            banner = line.split("tauren-md%0Av")[1].strip()
            break


def test_version_vs_readme():
    assert readme == version.__version__


def test_version_vs_banner():
    assert banner == version.__version__

    
def test_readme_vs_banner():
    assert readme == banner
    
