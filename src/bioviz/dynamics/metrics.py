import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import numpy as np
from typing import Dict, Any, Optional


class TrajectoryMetricsEngine:
    """
    Computes essential quantitative structural metrics across MD simulation trajectories.
    """

    def __init__(self, universe: mda.Universe):
        self.universe = universe

    # Global Root Mean Square Deviation (RMSD)
    def compute_rmsd(self, select_str: str = "protein and backbone", ref_frame: int = 0) -> Dict[str, np.ndarray]:
        print("[BioViz Metrics] Calculating global RMSD...")
        R = rms.RMSD(self.universe, self.universe, select=select_str, ref_frame=ref_frame)
        R.run()

        return {
            "frames": R.results.rmsd[:, 0],
            "time_ps": R.results.rmsd[:, 1],
            "rmsd": R.results.rmsd[:, 2]
        }

    # Root Mean Square Fluctuation (RMSF)
    def compute_rmsf(self, select_str: str = "protein and name CA") -> Dict[str, np.ndarray]:
        print("[BioViz Metrics] Calculating residual RMSF...")
        atom_group = self.universe.select_atoms(select_str)

        solver = rms.RMSF(atom_group)
        solver.run()

        return {
            "resnums": atom_group.resnums,
            "resnames": atom_group.resnames,
            "rmsf": solver.results.rmsf
        }
    # Dynamic hydrogen bond networks
    def compute_hydrogen_bonds(self, distance_cutoff: float = 3.0, angle_cutoff: float = 150.0) -> Dict[str, Any]:
        print("[BioViz Metrics] Analyzing hydrogen bond network dynamics...")

        # Initialize using accepted MDAnalysis signature parameters
        hb = HydrogenBondAnalysis(
            universe=self.universe,
            donors_sel="protein and (name N O)",
            hydrogens_sel="protein and name H*",
            acceptors_sel="protein and (name O or name N)",
            d_a_cutoff=distance_cutoff,
            d_h_a_angle_cutoff=angle_cutoff
        )
        hb.run()

        # Summary metrics per frame
        frame_counts = [len(ts_bonds) for ts_bonds in hb.results.hbonds] if hasattr(hb.results, 'hbonds') else []

        return {
            "total_timesteps_analyzed": len(hb.results.hbonds) if hasattr(hb.results, 'hbonds') else 0,
            "hbond_counts_per_frame": frame_counts,
            "average_hbonds": np.mean(frame_counts) if frame_counts else 0.0
        }