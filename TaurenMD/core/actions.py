actions_dict = {
    "remove_solvent": lambda x: x.remove_solvent(),
    "undo_remove_solvent": lambda x: x.undo_remove_solvent(),
    "align_traj": lambda x: x.align_traj(),
    "image_molecules": lambda x: x.image_molecules(),
    "select_frames": lambda x, y: x.select_frames(**y),
    "select_atoms": lambda x: x.select_atoms(**y),
    
    
    
    
    "frame_slice": lambda x, y: x.frame_slice(**y),
    "atom_selection": lambda x, y: x.set_atom_selection(**y),
    "align_traj": lambda x, y: x.align_traj(**y),
    "try_image_molecules": lambda x, y: x.image_molecules(**y),
    "frames2file": lambda x, y: x.frames2file(**y),
    "save_traj": lambda x, y: x.save_traj(**y),
    "produce_rmsds_combined_chains":
        lambda x, y: produce.rmsds_combined_chains(x, **y),
    "produce_rmsds_separated_chains":
        lambda x, y: produce.rmsds_separated_chains(x, **y),
    }
