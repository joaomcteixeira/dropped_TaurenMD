import string
import numpy as np

import MDAnalysis as mda
from MDAnalysis.analysis import align as mdaalign
from MDAnalysis.analysis.rms import RMSD as mdaRMSD

from TaurenMD import log
from TaurenMD.logger import (
    Title as T,
    Content as C,
    End as E,
    )

import TaurenMD.taurens.commons as TC

# from TaurenMD.taurens.abc_tauren import TaurenTraj

class TaurenMDAnalysis:  #(TaurenTraj):
    
    chain_or = ' or '
    
    def __init__(self, trajpath, topopath):
        self._universe = mda.Universe(topopath, trajpath)
        self._topology = mda.Universe(topopath)
        self.atom_selection = "all"
        self.chain_selection = 'all'
        self._total_n_frames = self.universe.trajectory.n_frames
        self._complete_frame_list = list(range(self._total_n_frames))
        self._frame_slicer = None
        self.solvent_selection = True
    
    @property
    def universe(self):
        return self._universe
    
    @property
    def topology(self):
        return self._topology
    
    @property
    def trajectory(self):
        return self.universe.select_atoms(self.selection)
    
    @property
    def total_n_frames(self):
        return self._total_n_frames
    
    @property
    def complete_frame_list(self):
        return self._complete_frame_list
    
    def select_frames(self, selection=None, start=None, end=None, step=None):
        
        
        if TC.is_selection_string_list(selection):
            a = [int(n) for n in selection.split()]
        elif TC.is_selection_list(selection):
            a = selection
        elif TC.is_selection_special_string(selection):
            a = TC.get_slice_from_special_string(selection)
        else:
            a = slice(start, end, step)
        
        self.frame_slicer = a
        return
    
    @property
    def frame_slicer(self):
        return self._frame_slicer
    
    @frame_slicer.setter
    def frame_slicer(self, slicer):
        valid = [
            slicer is None,
            isinstance(slicer, list),
            isinstance(slicer, slice),
            ]
        
        if any(valid):
            self._frame_slicer = slicer
        else:
            raise TypeError(f'slicer of wrong type, {type(slicer)}')
        
    
    @property
    def frames(self):
        """
        The frame list according to the current frame selection,
        :attr:`frame_slicer`.
        """
        
        if isinstance(self.frame_slicer, slice):
            stuple = self.frame_slicer.indices(len(self.complete_frame_list))
            list_of_frames = self.complete_frame_list[slice(*stuple)]
        
        elif self.frame_slicer is None  :
            list_of_frames = self.complete_frame_list
        
        elif isinstance(self.frame_slicer, list):
            list_of_frames = \
                [self.complete_frame_list[i] for i in self.frame_slicer]
        
        else:
            raise ValueError('could not get frames list')
        
        assert isinstance(list_of_frames, list), 'IS NOT A LIST!'
        return list_of_frames
    
    
    @property
    def n_frames(self):
        return len(self.frames)
    
    @property
    def pdb_name_fmt(self):
        return '{:0>' + str(len(str(len(self.frames)))) + '}'
    
    @property
    def atom_selection(self):
        """
        String that defines the current atom selection.
        
        If not defined assumes "all".
        """
        return self._atom_selection
    
    @atom_selection.setter
    def atom_selection(self, selector):
        """(str)"""
        
        if selector in (None, "all"):
            self._atom_selection = "all"
        
        elif not isinstance(selector, str):
            raise TypeError(
                "<selector> must be STRING type."
                f" '{type(selector)}' given."
                )
        
        else:
            self._atom_selection = f"({selector})"
    
    @property
    def chain_selection(self):
        return self._chain_selection
    
    @chain_selection.setter
    def chain_selection(self, chains):
        if chains in (None, 'all'):
            self._chain_selection = 'all'
            return
        
        elif isinstance(chains, str) \
                and all(c.isalnum() for c in chains.split(',')):
            chain_list = chains.split(',')
            
        elif isinstance(chains, list) and \
                all(isinstance(c, str) and c.isalnum() for c in chains):
            chain_list = chains
        
        chain_list = [f'segid {c}' for c in chain_list]
        
        self._chain_selection = self.chain_or.join(chain_list)
    
    @property
    def chain_list(self):
        
        if self.chain_selection == 'all':
            return [f'segid {x}' for x in self.filter_selectors(string.ascii_letters)]
        else:
            return self.chain_selection.split(self.chain_or)
    
    def filter_selectors(self, selectors_list):
        # https://www.mdanalysis.org/docs/documentation_pages/selections.html#simple-selections
        l = list(filter(
            lambda x: len(self.topology.select_atoms(f'segid {x}')) > 0,
            selectors_list,
            ))
        
        return l
    
    
    
    
    
    @property
    def solvent_selection(self):
        return self._solvent_selection
    
    @solvent_selection.setter
    def solvent_selection(self, selector):
        if selector is False:
            self._solvent_selection = "(protein or nucleic)"
        elif selector is True:
            self._solvent_selection = "all"
        elif isinstance(selector, list):
            self._solvent_selector = (
                "(protein or nucleic or "
                f"{'name' + ' or name '.join(selector)})"
                )
        else:
            raise ValueError("solvent selection could not be assigned.")
    
    @property
    def selection(self):
        
        sel = "({}) and ({}) and ({})".format(
            self.atom_selection,
            self.chain_selection,
            self.solvent_selection,
            )
        
        return sel
    
    def remove_solvent(self, include=False, **kwargs):
        self.solvent_selection = include
        return
    
    def undo_rmv_solvent(self):
        self.solvent_selection = True
        return
    
    def frames2files(self, prefix='_', ext='pdb'):
        for frame in self.frames:
            
            file_name = '{}{}.{}'.format(
                prefix,
                self.pdb_name_fmt.format(frame),
                ext,
                )
            
            self.trajectory.write(
                filename=file_name,
                frames=[frame],
                )
            log.info(C(f'saved {file_name}'))
        return
    
    def frames2models(self, filename='mytraj.pdb'):
        self.trajectory.write(
            filename=filename,
            frames=self.frames,
            )
        return
    
    def save_traj(self, filename='mytraj.dcd'):
        selection = self.universe.select_atoms(self.selection)
        with mda.Writer(filename, selection.n_atoms) as W:
            for ts in self.trajectory[self.frames]:
                W.write(selection)
                log.info(C(f"exported {ts}"))
    
    def align_traj(self, weight='mass'):
        
        alignment = mdaalign.AlignTraj(
            self.universe,
            self.topology,
            in_memory=True,
            verbose=True,
            weight=weight,
            )
        
        alignment.run()
        
        return
    
    def calc_rmsds_combined_chains(self):
        
        log.info(C('for selection: {}'.format(self.selection)))
        log.info(C('for frames: {}'.format(self.frame_slicer)))
        
        R = mdaRMSD(
            self.universe,
            self.topology,
            select=self.selection,
            groupselection=None,
            ref_frame=0,
            verbose=False,
            )
        
        R.run(verbose=False)
        
        # rmsds[:, ii] = R.rmsd[:, 2][self._fslicer]
        
        return R.rmsd[self.frames, 2]
    
    def calc_rmsds_separated_chains(self):
        
        rmsds = np.empty((self.n_frames, len(self.chain_list)))
        
        for ii, chain in enumerate(self.chain_list):
            
            final_selection = (
                f"{self.atom_selection}"
                f" and ({chain})"
                f" and ({self.solvent_selection})"
                )
            
            log.info(C(f'calculating RMSDs for {final_selection}'))
            
            R = mdaRMSD(
                self.universe.select_atoms(final_selection),
                self.topology.select_atoms(final_selection),
                groupselection=None,
                ref_frame=0,
                verbose=False,
                )
            
            R.run(verbose=False)
            
            rmsds[:, ii] = R.rmsd[self.frames, 2]
        
        return rmsds.T
        
            # subplot_has_data.append(True)
        
        # column_headers = list(map(
            # lambda x: x.replace("segid ", ""),
            # filtered_selectors
            # ))
        
        # assert isinstance(column_headers, list), "c_selectors NOT list!"
        
        # return (
            # rmsds[:, subplot_has_data],
            # np.array(column_headers)[subplot_has_data],
            # )
