import os   

from tauren import produce, load
from tests import commons

file_path = "tests"
rf = "reference"
mdtype = "mda"
trajtype = "mdanalysis"

trajectory = os.path.join(file_path, rf, "traj_test_PCNA.dcd")
topology = os.path.join(file_path, rf, "topology_test.pdb")

rmsd_noHOH = os.path.join(
    file_path,
    rf,
    mdtype,
    "rmsds_combined_chains_-protein_or_nucleic-_all_all.csv"
    )

rmsd_sep_noHOH = os.path.join(
    file_path,
    rf,
    mdtype,
    "rmsds_separated_chains_-protein_or_nucleic-_all_all.csv"
    )

traj = load.load_traj(
    trajectory,
    topology,
    traj_type=trajtype,
    )

traj.remove_solvent()
traj.calc_rmsds_combined_chains()
traj.calc_rmsds_separated_chains()

def test_results_len():
    
    assert len(traj.observables) == 2

def test_export_combined_1():
    
    export_data = {
        "file_name": None,
        "sep": ","
        }
    
    produce._update_export_data(
        export_data,
        traj.observables[0],
        0,
        )
    
    assert export_data["file_name"] == \
        "data_rmsds_combined_chains_-protein_or_nucleic-_all_all"

def test_export_combined_2():
    
    export_data = {
        "file_name": "mydata.csv",
        "sep": ","
        }
    
    produce._update_export_data(
        export_data,
        traj.observables[0],
        0,
        )
    
    assert export_data["file_name"] == "mydata.csv"

def test_export_combined_3():
    
    export_data = {
        "file_name": None,
        "sep": ",",
        "suffix" : "csv",
        }
    
    produce._update_export_data(
        export_data,
        traj.observables[0],
        0,
        )
    
    assert export_data["file_name"] == \
        "data_rmsds_combined_chains_-protein_or_nucleic-_all_all"


def test_columns_1():
    
    assert traj.observables[0]["columns"] == ["frames", "A,B,C"]

def test_columns_2():
    
    assert traj.observables[1]["columns"] == ["frames", "A", "B" , "C"]
    
def test_update_single_plot_1():
    
    kw = {
        "label": "my label 1",
        "fig_name": None,
            }
    
    produce._update_single_plot_config(
        kw,
        0,
        "plot",
        traj.observables[0],
        )
    
    assert kw["label"] == "my label 1"
    assert kw["fig_name"] == \
        "plot_rmsds_combined_chains_-protein_or_nucleic-_all_all"

def test_update_single_plot_2():
    
    kw = {
        "label": None,
        "fig_name": "myfig.pdf",
            }
    
    produce._update_single_plot_config(
        kw,
        0,
        "plot",
        traj.observables[0],
        )
    
    assert kw["label"] == "A,B,C"
    assert kw["fig_name"] == "myfig.pdf"

def test_update_multiple_plot_1():
    
    kw = {
        "labels": [
            "my label 1",
            "my label 2",
            "my label 3",
            ],
        "fig_name": None,
        "colors": None,
            }
    
    produce._update_multiple_plot_config(
        kw,
        1,
        "plot",
        traj.observables[1],
        )
    
    assert kw["labels"] == ["my label 1", "my label 2", "my label 3"]
    assert kw["fig_name"] == \
        "plot_rmsds_separated_chains_-protein_or_nucleic-_all_all"
    assert "colors" not in kw

def test_update_multiple_plot_2():
    
    kw = {
        "labels": None,
        "fig_name": "myfig.pdf",
        "colors": ["k"],
            }
    
    produce._update_multiple_plot_config(
        kw,
        1,
        "plot",
        traj.observables[1],
        )
    
    assert kw["labels"] == ["A", "B", "C"]
    assert kw["fig_name"] == "myfig.pdf"
    assert "colors" in kw

def test_produce_combined_1():
    
    calc_rmsds_combined_chains = {
        "chains": "all",
        "ref_frame": 0
        }
        
    export_data = {
        "file_name": None,
        "sep": ",",
        "suffix": "csv",
        "tojson": False
        }
        
    plot_rmsd_combined_chains  = {
        "label": None,
        "suptitle": "Combined Chain RMSDs",
        "x_label": "Frame Number",
        "y_label": "RMSDs",
        "color": "blue",
        "alpha": 0.7,
        "grid": True,
        "grid_color": "lightgrey",
        "grid_ls": "-",
        "grid_lw": 1,
        "grid_alpha": 0.5,
        "legend": True,
        "legend_fs": 6,
        "legend_loc": 4,
        "fig_name": None
        }
    
    produce.rmsds_combined_chains(
        traj,
        calc_rmsds_combined_chains,
        export_data=export_data,
        plot_rmsd_combined_chains=plot_rmsd_combined_chains,
        )
    
    rmsdtest = "data_rmsds_combined_chains_-protein_or_nucleic-_all_all.csv"
    
    assert commons.compare_csv(rmsd_noHOH, rmsdtest)
    os.remove(rmsdtest)
    os.remove("plot_rmsds_combined_chains_-protein_or_nucleic-_all_all.pdf")

def test_produce_combined_1():
    
    calc_rmsds_separated_chains = {
        "chains": "all",
        "ref_frame": 0
        }
        
    export_data = {
        "file_name": None,
        "sep": ",",
        "suffix": "csv",
        "tojson": False
        }
        
    plot_rmsd_chain_per_subplot  = {
        "labels": None,
        "suptitle": "RMSDs per chain",
        "x_label": "Frame Number",
        "y_label": "RMSDs",
        "colors": None,
        "alpha": 0.7,
        "grid": True,
        "grid_color": "lightgrey",
        "grid_ls": "-",
        "grid_lw": 1,
        "grid_alpha": 0.5,
        "legend": True,
        "legend_fs": 6,
        "legend_loc": 4,
        "fig_name": None
        }
    
    plot_rmsd_individual_chains_one_subplot = {
        "labels": None,
        "suptitle": "Chains' RMSDs",
        "x_label": "Frame Number",
        "y_label": "RMSDs",
        "colors": None,
        "alpha": 0.7,
        "grid": True,
        "grid_color": "lightgrey",
        "grid_ls": "-",
        "grid_lw": 1,
        "grid_alpha": 0.5,
        "legend": True,
        "legend_fs": 6,
        "legend_loc": 4,
        "fig_name": None
        }
    
    produce.rmsds_separated_chains(
        traj,
        calc_rmsds_separated_chains,
        export_data=export_data,
        plot_rmsd_chain_per_subplot=plot_rmsd_chain_per_subplot,
        plot_rmsd_individual_chains_one_subplot=plot_rmsd_individual_chains_one_subplot,
        )
    
    rmsdtest = "data_rmsds_separated_chains_-protein_or_nucleic-_all_all.csv"
    
    assert commons.compare_csv(rmsd_sep_noHOH, rmsdtest)
    os.remove(rmsdtest)
    os.remove("plot_rmsd_chain_per_subplot_rmsds_separated_chains_-protein_or_nucleic-_all_all.pdf")
    os.remove("plot_rmsd_individual_chains_one_subplot_rmsds_separated_chains_-protein_or_nucleic-_all_all.pdf")
        
