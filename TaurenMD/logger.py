import sys
import logging
import logging.handlers

log_filename = 'taurenmd.log'
debug_filename = 'taurenmd.debug'


class Title:
    
    def __init__(self, msg):
        self.m = msg
    
    def __str__(self):
        return f"* {self.m.title()} ..."

class Content:
    
    def __init__(self, msg):
        self.m = msg
    
    def __str__(self):
        return f"    - {self.m}"


class End:
    
    def __init__(self, msg):
        self.m = msg
    
    def __str__(self):
        return f"    ...{self.m}"


def get_logger(name):
    """Returns a configured logger."""
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    if not logger.handlers:
        # debugger
        dbgfmt = logging.Formatter(
            "%(message)s"
            "`%(levelname)s - "
            "%(filename)s:%(name)s:%(funcName)s:%(lineno)d`\n"
            )
        debughandler = logging.handlers.RotatingFileHandler(
            "tauren.debug",
            maxBytes=10485760,
            mode="a",
            )
        debughandler.setFormatter(dbgfmt)
        debughandler.setLevel(logging.DEBUG)
        
        debugmemhandler = logging.handlers.MemoryHandler(
            1000,
            target=debughandler,
            flushLevel=logging.ERROR,
            )
        
        debugmemhandler.setFormatter(dbgfmt)
        debugmemhandler.setLevel(logging.DEBUG)
        
        # user info log
        nfofmt = logging.Formatter("%(message)s")
        infohandler = logging.handlers.RotatingFileHandler(
            "tauren.log",
            maxBytes=10485760,
            mode="a",
            )
        
        infohandler.setLevel(logging.INFO)
        infohandler.setFormatter(nfofmt)
        
        infomemhandler = logging.handlers.MemoryHandler(
            1000,
            target=infohandler,
            flushLevel=logging.ERROR,
            )
        
        infomemhandler.setLevel(logging.INFO)
        infomemhandler.setFormatter(nfofmt)
        
        # console
        ch = logging.StreamHandler(stream=sys.stdout)
        ch.setLevel(logging.INFO)
        ch.setFormatter(nfofmt)
        
        # add handlers
        logger.addHandler(debugmemhandler)
        logger.addHandler(infomemhandler)
        logger.addHandler(ch)
    
    return logger
