import os
import contextlib

from TaurenMD.logger import debug_filename, log_filename
from TaurenMD.logger import get_logger

with contextlib.suppress(FileNotFoundError):
    os.remove(debug_filename)
    os.remove(log_filename)

log =  get_logger(__name__)

from TaurenMD.base import TaurenMD

# from TaurenMD.taurens.taurenmda import TaurenMDAnalysis
