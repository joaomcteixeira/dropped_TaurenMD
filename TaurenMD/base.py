from pathlib import Path

from TaurenMD import log
from TaurenMD.logger import (
    Title as T,
    Content as C,
    End as E,
    )

from TaurenMD.taurens import trajtypes
import TaurenMD.plots as TPLOT

class TaurenMD:
    """
    TaurenMD interface.
    """
    
    def __init__(self, trajpath, topopath, mdlib='mdanalysis'):
        
        self.reload_system(trajpath, topopath, mdlib)
        
        return
    
    def __str__(self):
        return repr(self)
    
    def __repr__(self):
        
        r = (
            f"TaurenMD(trajpath='{self.trajpath.resolve()}',"
            f"topopath='{self.topopath.resolve()}',"
            f"mdlib='{self.mdlib}')"
            )
        return r
    
    def __bool__(self):
        return bool(self.taurentraj)
    
    def reload_system(self, trajpath, topopath, mdlib):
        
        self.trajpath = Path(trajpath)
        self.topopath = Path(topopath)
        self.mdlib = mdlib
        self.taurentraj = trajtypes[mdlib](trajpath, topopath)
        return
    
    def remove_solvent(self, **kwargs):
        """
        Removes solvent from trajectory.
        
        Parameters
        ----------
        include : :obj:`list` or :obj:`bool`
            Available only when using :mdtrajdoc:`MDTraj <>`
            as :doctrajtype:`trajectory type <>`.
            List of solvent residue names to retain in the new
            trajectory.
            Defaults to False, no solvent residues are kept.
            :mdtdocrmvsol:`MDTraj remove_solvent documentation <>`.
        """
        log.info(T("removing solvent"))
        self.taurentraj.remove_solvent(**kwargs)
        log.info(E("done"))
    
    def undo_remove_solvent(self):
        """
        Undo a previous action of solvent removal by activating
        solvent atoms again.
        """
        self.taurentraj.undo_rmv_solvent()
    
    def align_traj(self):
        """
        Aligns the whole trajectory to its topology.
        
        Currently only implemented for :mdanalysis:`MDAnalysis <>`
        subroutines.
        """
        log.info(T("aligning trajectory"))
        self.taurentraj.align_traj()
    
    def image_molecules(self, **kwargs):
        """
        Images molecules.
        
        Currently only available for :mdtrajdoc:`MDTraj <>`
        subroutines.
        
        Images molecules according to :mdjdocim:`MDTraj.image_molecules <>`.
        Receives the same arguments.
        """
        log.info(T("imaging molecules"))
        self.taurentraj.image_molecules(**kwargs)
    
    def select_frames(self, selection=None, start=None, end=None, step=None):
        """
        Slices trajectory in frames for the subsequent operations.
        
        Parameters
        ----------
        selection: :obj:`str`
            Special selection string for trajectory slicing.
            Acceptted formats:
            - "START"
            - "START:"
            - "START:END"
            - "START:END:STEP"
            - "::STEP"
        
        start : int
            The starting frame.
            Frame index starts at 0.
            Defaults to None, slices from the first frame.
        
        end : int
            The end frame for the new slicing (exclusive).
            Defaults to None, slices to the last frame.
        
        step : int
            Integer value which determines the increment between
            each index for slicing.
            Defaults to None, equals to 1.
        """
        log.info(
            T("slicing trajectory: [{}, {}, {}, {}]".format(
                selection,
                start,
                step,
                end,
                ),
            ),
            )
        self.taurentraj.select_frames(selection, start, end, step)
    
    def select_atoms(self, selector=None):
        """
        Sets the atom selection.
        
        Atom selection will be used in subsequent operations.
        
        Parameters
        ----------
        selector : :obj:`str`.
            The selection string.
            If None type provided assumes ``"all"``.
            Should be of a valid format according to the
            MD analysis library used in :doctrajtype:`trajectory_type <>`.
            Please refer to the MD analysis library
            specific documentation:
            
            - :mdaselections:`MDAnalysis <>`
            - :mdtselections:`MDTraj <>`
        """
        log.info(T("setting atom selection"))
        self.taurentraj.atom_selection = selector
        log.info(
            C("atom selection set to %s"),
            self.taurentraj.atom_selection,
            )
        return
    
    def select_chains(self, chains):
        """
        Selects chains.
        
        Parameter
        ---------
        chains : obj:`str` or :obj:`list`
            A string with the name of the chain or a list with chain names.
            Defaults to None, select all chains.
        """
        log.info(T('selecting chains'))
        self.taurentraj.chain_selection = chains
        log.info(E('done'))
    
    def frames2files(self, prefix="_", ext="pdb"):
        """
        Exports selected frames to files according to <ext>.
        
        To select frames use .select_frames method.
        """
        log.info(T("extracting frames"))
        self.taurentraj.frames2files(prefix, ext)
        log.info(E('done'))
        return
    
    def frames2models(self, filename='mytraj.pdb'):
        """
        Exports selected frames to models in a single file.
        
        File type is deduced from filename extension.
        
        To select frames use .select_frames method.
        """
        log.info(T('extracting frames to models'))
        log.info(C(filename))
        self.taurentraj.frames2models(filename)
        log.info(E('done'))
        return
    
    def save_traj(self, file_name="traj_output.dcd", **kwargs):
        """
        Saves trajectory to file.
        
        ATTENTION: Overwrites existing files.
    
        Parameters
        ----------
        file_name : :obj:`str`
            Name of the output trajectory file.
            File extension is taken from ``file_name``.
        """
        log.info(T("exporting trajectory to"))
        log.info(C(file_name))
        self.taurentraj.save_traj(file_name)
        log.info(E("saved"))
        return

    def calc_rmsds_combined_chains(
            self,
            plot=False,
            plot_file_name='RMSDs_for_chains_combined.pdf',
            ):
        """
        Calculates the global RMSD by considering all chains combined.
        """
        log.info(T("calculating rmsd for combined chains"))
        self.latest_results = self.taurentraj.calc_rmsds_combined_chains()
        log.info(E("done"))
        
        if plot:
            log.info(T('plotting parameter'))
            TPLOT.PP(
                self.taurentraj.frames,
                self.latest_results,
                fig_name=plot_file_name,
                **TPLOT.DEFS.rmsd_combined_chains,
                )
            log.info(E(f'plotted {plot_file_name}'))
        return
    
    def calc_rmsds_separated_chains(
            self,
            plot_single_plot=False,
            pss_fig_name='RMSDs_separated_chains_single_plot.pdf',
            plot_param_series_subplots=False,
            ):
        
        log.info(T("calculating rmsd for chains individually"))
        self.latest_results = self.taurentraj.calc_rmsds_separated_chains()
        log.info(E("RMSDs calculated"))
        
        if plot_single_plot:
            log.info(T('plotting series in multiple subplots'))
            TPLOT.PPSS(
                self.taurentraj.frames,
                self.latest_results,
                fig_name=filename,
                **TPLOT.DEFS.rmsd_chain_per_subplot
                )
            log.info(E(f'plotted {filenme}'))
            # plot_series_single_subplot(filename='RMSDs_separated_chains_single_plot.pdf')
        
        if plot_param_series_subplots:
            # plot_series_subplots(filename='RMSDs_separated_chains.pdf')
            log.info(T('plotting series in single subplot'))
            TPLOT.PPSSS(
                self.taurentraj.frames,
                self.latest_results,
                fig_name=filename,
                **TPLOT.DEFS.rmsd_individual_chains_one_subplot,
                )
            
            log.info(E(f'plotted {filename}'))
        return
    
    # def plot_RMSD_combined_chains(self, filename='plot_param.pdf'):
        
        # log.info(T('plotting parameter'))
        # TPLOT.PP(
            # self.taurentraj.frames,
            # self.latest_results,
            # fig_name=file_name,
            # **TPLOT.DEFS.rmsd_combined_chains,
            # )
        # log.info(E(f'plotted {filename}'))
    
    # def plot_series_subplots(self, filename='series_subplots.pdf'):
        
        
    
    # def plot_series_single_subplot(self, filename='series_single_plot.pdf'):
        
        # log.info(T('plotting series in single subplot'))
        # TPLOT.PPSSS(
            # self.taurentraj.frames,
            # self.latest_results,
            # fig_name=filename,
            # **TPLOT.DEFS.rmsd_individual_chains_one_subplot,
            # )
        
        # log.info(E(f'plotted {filename}'))
