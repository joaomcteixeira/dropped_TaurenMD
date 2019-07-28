import itertools as it

from matplotlib import pyplot as plt
from matplotlib import colors as mcolors

from TaurenMD import log
from TaurenMD.logger import (
    Title as T,
    Content as C,
    End as E,
    )

import TaurenMD.plots.commons as TPC

def plot_param_series_single_subplot(
        x_data,
        y_data,
        *,
        labels=None,
        suptitle=None,
        x_label=None,
        y_label=None,
        colors=list(mcolors.BASE_COLORS.keys()),
        alpha=0.7,
        grid=True,
        grid_color="lightgrey",
        grid_ls="-",
        grid_lw=1,
        grid_alpha=0.5,
        legend=True,
        legend_fs=6,
        legend_loc=4,
        fig_name="plot_param_series_single_subplot.pdf",
        **kwargs
        ):
    """
    Plots a single plot with the combined RMSD.
    
    Bellow parameters concern data representation and are considered
    of highest importance because their incorrect use can mislead
    data analysis and consequent conclusions.
    
    Plot style parameters concernning only plot style, i.e., colors,
    shapes, fonts, etc... and which do not distort the actual data,
    are not listed in the paremeter list bellow. We hope these
    parameter names are self-explanatory and are listed in the function
    definition.
    
    Parameters
    ----------
    x_data : interable
        Contains X axis data.
        Should be accepted by matplotlib.
    
    y_data : np.ndarray, shape=(N, M)
        Container of the Y axis data.
        Where N is the number of chains (data series), and
        M the data for each series (RMSDs of a given chain)
    
    labels : str, optional
        The chain labels to represent in plot legend.
        Defauts to: "no labels provided".
    
    fig_name : str, optional
        The file name with which the plot figure will be saved
        in disk. Defaults to rmsd_individual_chains_one_subplot.pdf.
        You can change the file type by specifying its extention in
        the file name.
    """
    log.info("* Plotting chains RMSDs single subplot...")
    
    fig, ax = plt.subplots(nrows=1, ncols=1)
    
    plt.tight_layout(rect=[0.05, 0.02, 0.995, 0.985])
    fig.suptitle(
        suptitle,
        x=0.5,
        y=0.990,
        va="top",
        ha="center",
        )
    
    colors_it = it.cycle(colors)
    max_rmsds = 0
    
    if not labels:
        labels = list(range(y_data.shape[0]))
    
    for ii, chain_rmsd_data in enumerate(y_data):
        
        ax.plot(
            x_data,
            chain_rmsd_data,
            label=labels[ii],
            color=next(colors_it),
            alpha=alpha,
            )
        
        if chain_rmsd_data.max() > max_rmsds:
            max_rmsds = chain_rmsd_data.max()
    
    ax.set_xlim(x_data[0], x_data[-1])
    ax.set_ylim(0, max_rmsds * 1.1)
    
    ax.set_xlabel(x_label, weight='bold')
    ax.set_ylabel(y_label, weight='bold')
    
    if grid:
        ax.grid(
            color=grid_color,
            linestyle=grid_ls,
            linewidth=grid_lw,
            alpha=grid_alpha,
            )
    
    if legend:
        ax.legend(
            fontsize=legend_fs,
            loc=legend_loc,
            )
    
    fig.savefig(fig_name)
    log.info(TPC.msg_fig_saved.format(fig_name))
    
    plt.close("all")
    
    return