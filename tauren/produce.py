"""
Functions that perform combined operations that are related.
"""
# Copyright © 2018-2019 Tauren-MD Project
#
# Tauren-MD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Tauren-MD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Tauren-MD. If not, see <http://www.gnu.org/licenses/>.
#
# Contributors to this file:
# - João M.C. Teixeira (https://github.com/joaomcteixeira)
from tauren import logger, core, plot

log = logger.get_log(__name__)


@core.log_args
def rmsds_combined_chains(
        taurentraj,
        calc_rmsds_combined_chains,
        *,
        export_data=False,
        plot_rmsd_combined_chains=False,
        **kwargs
        ):
    """
    Execute routines related to RMSDs of combined chains.
    """
        
    index = taurentraj.calc_rmsds_combined_chains(
        **calc_rmsds_combined_chains
        )
    
    if export_data:
        
        _update_export_data(
            export_data,
            taurentraj.observables[index],
            index,
            )
        
        log.debug(export_data)
        
        taurentraj.export_data(index, **export_data)
    
    if plot_rmsd_combined_chains:
        
        _update_single_plot_config(
            plot_rmsd_combined_chains,
            index,
            "plot",
            taurentraj.observables[index],
            )
        
        log.debug(plot_rmsd_combined_chains)
        
        plot.rmsd_combined_chains(
            taurentraj.observables[index]["data"][:, 0],
            taurentraj.observables[index]["data"][:, 1],
            **plot_rmsd_combined_chains,
            )
        
    return


@core.log_args
def rmsds_separated_chains(
        taurentraj,
        calc_rmsds_separated_chains,
        *,
        export_data=False,
        plot_rmsd_chain_per_subplot=False,
        plot_rmsd_individual_chains_one_subplot=False,
        **kwargs
        ):
    """
    Execute routines related to RMSDs of separated chains.
    """
        
    index = taurentraj.calc_rmsds_separated_chains(
        **calc_rmsds_separated_chains
        )
    
    if export_data:
        
        _update_export_data(
            export_data,
            taurentraj.observables[index],
            index,
            )
        
        log.debug(export_data)
        
        taurentraj.export_data(index, **export_data)
    
    if plot_rmsd_chain_per_subplot:
        
        _update_multiple_plot_config(
            plot_rmsd_chain_per_subplot,
            index,
            "plot_rmsd_chain_per_subplot",
            taurentraj.observables[index],
            )
        
        log.debug(plot_rmsd_chain_per_subplot)
        
        plot.rmsd_chain_per_subplot(
            taurentraj.observables[index]["data"][:, 0],
            taurentraj.observables[index]["data"][:, 1:].T,
            **plot_rmsd_chain_per_subplot,
            )
    
    if plot_rmsd_individual_chains_one_subplot:
        
        _update_multiple_plot_config(
            plot_rmsd_individual_chains_one_subplot,
            index,
            "plot_rmsd_individual_chains_one_subplot",
            taurentraj.observables[index],
            )
        
        log.debug(plot_rmsd_individual_chains_one_subplot)
        
        plot.rmsd_individual_chains_one_subplot(
            taurentraj.observables[index]["data"][:, 0],
            taurentraj.observables[index]["data"][:, 1:].T,
            **plot_rmsd_individual_chains_one_subplot,
            )
        
    return


def _get_key_list(key):
    """
    .. deprecated:: 0.6.0
    """
    kl = key.split(",")
    log.debug(kl)
    return kl


@core.log_args
def _update_export_data(
        kwargs,
        data_dict,
        index,
        ):
        
    if kwargs["file_name"] is None:
        try:
            kwargs["file_name"] = f"data_{data_dict['name']}"
        
        except KeyError:
            kwargs["file_name"] = f"data_for_index_{index}"


@core.log_args
def _update_single_plot_config(
        kwargs,
        index,
        name,
        data_dict,
        ):
        
    if kwargs["label"] is None:
        kwargs["label"] = ", ".join(data_dict["columns"][1:])
    
    if kwargs["fig_name"] is None:
        try:
            kwargs["fig_name"] = f"{name}_{data_dict['name']}"
        
        except KeyError:
            kwargs["fig_name"] = f"{name}_data_for_index_{index}"
        

@core.log_args
def _update_multiple_plot_config(
        kwargs,
        index,
        name,
        data_dict,
        ):
    
    if kwargs["labels"] is None:
        kwargs["labels"] = data_dict["columns"][1:]
    
    if kwargs["fig_name"] is None:
        try:
            kwargs["fig_name"] = f"{name}_{data_dict['name']}"
        
        except KeyError:
            kwargs["fig_name"] = f"{name}_data_for_index_{index}"
    
    if kwargs["colors"] is None:
        kwargs.pop("colors")
