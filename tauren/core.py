"""
Variables and functions that serve system-wide.
"""
import functools
from tauren import logger

log = logger.get_log(__name__)

trajectory_types = (".xtc", ".nc", ".trr", ".h5", ".pdb", ".binpos", ".dcd")
"""Types os trajectories accepted."""

topology_types = (".pdb", ".cif")
"""Types of topologies accepted."""


def log_args(func):
    
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        
        log.debug(f"LOG ARGUMENTS FOR FUNC: {func.__name__}")
        
        for ii, arg in enumerate(args, start=1):
            log.debug(f"\tpositional arg {ii}: {repr(arg)} of type {type(arg)}")
        
        for key, value in kwargs.items():
            
            log.debug(f"\tkw {key}: {repr(value)} of type {type(value)}")
        
        result = func(*args, **kwargs)
        
        log.debug(f"\tRETURNS: {result} of type {type(result)}")
        
        return result
    
    return wrapper
