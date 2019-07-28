from TaurenMD import log

msg_fig_saved = "* Saved figure: {}"

def calc_fig_size(
        nsubplots,
        *,
        ncols=1,
        irow=4.8,
        icol=6.4,
        ):
    """
    Calculates the figure dimensions (size in inches).
    
    Figure dimensions are calculated based on the following parameters:
    
    Parameters
    ----------
    nsubplots : int
        The total number of subplots
    
    ncols : int, optional, defaults 1
        The desired number of columns
    
    irow : float, optional, defaults 4.8
        Number os inches per row.
    
    icol : float, optional, defaults 6.4
        Number of inches per column.
    
    
    Returns
    -------
    A tuple (width, height) in inches.
    """
    
    # https://matplotlib.org/api/_as_gen/matplotlib.pyplot.figure.html
    
    width = ncols * icol
    height = (nsubplots // ncols) * irow
    
    fs = (width, height)
    
    log.debug(f"returning: {fs}")
    
    return fs