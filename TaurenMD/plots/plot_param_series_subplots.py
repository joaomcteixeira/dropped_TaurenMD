import itertools as it
import numpy as np

from matplotlib import pyplot as plt
from matplotlib import colors as mcolors

from TaurenMD import log
from TaurenMD.logger import (
    Title as T,
    Content as C,
    End as E,
    )

import TaurenMD.plots.commons as TPC 

def plot_param_series_subplots(
        x_data,
        y_data,
        *,
        labels="No labels provided",
        suptitle="RMSDs per chain",
        x_label="Frame Number",
        y_label="RMSDs",
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
        fig_name="rmsd_chain_per_subplot.pdf",
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
        Container of the X axis data. Should be accepted
        by matplotlib.
    
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
    log.info("* Plotting RMSDs per chain...")
    
    figsize = TPC.calc_fig_size(
        y_data.shape[0],
        ncols=1,
        irow=3,
        )
    
    fig, axs = plt.subplots(
        nrows=y_data.shape[0],
        ncols=1,
        figsize=figsize,
        sharex=True,
        )
    
    fig.suptitle(
        suptitle,
        x=0.5,
        y=0.990,
        va="top",
        ha="center",
        )
    
    try:
        axs = axs.ravel()
    
    except AttributeError:  #in case there only one subplot
        axs = np.array([axs])
    
    plt.tight_layout(rect=[0.05, 0.02, 0.995, 0.985])
    fig.subplots_adjust(hspace=0)
    
    colors_it = it.cycle(colors)
    
    max_rmsd = 0
    for ii, chain_rmsds in enumerate(y_data):
        
        axs[ii].plot(
            x_data,
            chain_rmsds,
            label=labels[ii],
            color=next(colors_it),
            alpha=alpha,
            )
        
        axs[ii].set_xlim(x_data[0], x_data[-1])
        
        axs[ii].set_ylabel(y_label, weight='bold')
        
        if grid:
            axs[ii].grid(
                color=grid_color,
                linestyle=grid_ls,
                linewidth=grid_lw,
                alpha=grid_alpha,
                )
    
        if legend:
            axs[ii].legend(
                fontsize=legend_fs,
                loc=legend_loc,
                )
        
        if chain_rmsds.max() > max_rmsd:
            max_rmsd = chain_rmsds.max()
        
    else:
        axs[ii].set_xlabel(x_label, weight='bold')
    
    log.debug(f"<max_rmsds>: {max_rmsd}")
    
    for axis in axs:
        
        axis.set_ylim(0, max_rmsd * 1.1)
        all_ticks = axis.get_yticks()
        axis.set_yticks(all_ticks[1:-1])
    
    else:
        axis.set_yticks(all_ticks[:-1])
    
    fig.savefig(fig_name)
    log.info(TPC.msg_fig_saved.format(fig_name))
    
    plt.close("all")
    
    return