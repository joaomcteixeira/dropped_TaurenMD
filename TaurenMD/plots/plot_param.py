from matplotlib import pyplot as plt

from TaurenMD import log
from TaurenMD.logger import (
    Title as T,
    Content as C,
    End as E,
    )
from TaurenMD.plots.commons import msg_fig_saved

# from TaurenMD.plots import commons as PLTC

def plot_param(
        x_data,
        y_data,
        *,
        label="No labels provided",
        suptitle=None,
        x_label=None,
        y_label=None,
        color='blue',
        alpha=0.7,
        grid=True,
        grid_color="lightgrey",
        grid_ls="-",
        grid_lw=1,
        grid_alpha=0.5,
        legend=True,
        legend_fs=6,
        legend_loc=4,
        fig_name='plot_param.pdf',
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
    x_data : interable of numbers
        Container of the X axis data. Should be accepted
        by matplotlib.
    
    y_data : np.ndarray, shape=(M,)
        Container of the Y axis data.
        Where M is the RMSDs data for the combined chains.
    
    label : str, optional
        The label to represent in plot legend.
        Defauts to: "no labels provided".
    
    fig_name : str, optional
        The file name with which the plot figure will be saved
        in disk. Defaults to rmsd_individual_chains_one_subplot.pdf.
        You can change the file type by specifying its extention in
        the file name.
    """
    # log.info("* Plotting combined Chain RMSDs...")
    
    fig, ax = plt.subplots(nrows=1, ncols=1)
    
    plt.tight_layout(rect=[0.05, 0.02, 0.995, 0.985])
    
    fig.suptitle(
        suptitle,
        x=0.5,
        y=0.990,
        va="top",
        ha="center",
        )
    
    ax.plot(
        x_data,
        y_data,
        label=label,
        color=color,
        alpha=alpha,
        )
    
    ax.set_xlabel(x_label, weight='bold')
    ax.set_ylabel(y_label, weight='bold')
    
    ax.set_xlim(x_data[0], x_data[-1])
    ax.set_ylim(0)
    
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
    log.info(C(msg_fig_saved.format(fig_name)))
        
    plt.close("all")
    
    return