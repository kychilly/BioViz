"""
BioViz Dynamics - Ingestion and Alignment Engine
Handles loading of atomic topologies and coordinate trajectories,
automatic PBC unwrapping, centering, and least-squares backbone alignment.
"""

import os
from typing import Optional
import MDAnalysis as mda
from MDAnalysis.transformations import unwrap, center_in_box
from MDAnalysis.analysis import align

class TrajectoryIngestionEngine:
    """
    Ingests molecular dynamics topologies and trajectories, standardizes
    periodic boundary conditions, centers systems, and aligns backbones.
    """

    def __init__(self, topology_path: str, trajectory_path: str):
        self.topology_path = topology_path
        self.trajectory_path = trajectory_path
        self.universe: Optional[mda.Universe] = None
        self._load_universe()

    # Loads topologies
    def _load_universe(self) -> None:
        if not os.path.exists(self.topology_path):
            raise FileNotFoundError(f"Topology file not found: {self.topology_path}")
        if not os.path.exists(self.trajectory_path):
            raise FileNotFoundError(f"Trajectory file not found: {self.trajectory_path}")

        try:
            self.universe = mda.Universe(self.topology_path, self.trajectory_path)
            print(f"[BioViz] Successfully loaded Universe: {len(self.universe.atoms)} atoms, {len(self.universe.trajectory)} frames.")
        except Exception as e:
            raise RuntimeError(f"Failed to load trajectory files into MDAnalysis Universe: {e}")

    def preprocess_trajectory(self, output_xtc_path: Optional[str] = None) -> mda.Universe:
        # MD cleanup
        if self.universe is None:
            raise ValueError("Universe is not initialized.")

        protein = self.universe.select_atoms("protein")
        if len(protein) == 0:
            target_selection = self.universe.select_atoms("all")
        else:
            target_selection = protein

        print("[BioViz] Running preprocessing pipeline...")

        # PBC transformations
        try:
            if self.universe.trajectory.ts.dimensions is not None:
                transformations = [
                    unwrap(target_selection),
                    center_in_box(target_selection, wrap=True)
                ]
                self.universe.trajectory.add_transformations(*transformations)
                print("[BioViz] PBC unwrap & centering applied successfully.")
            else:
                print("[BioViz Notice] No periodic box dimensions detected in trajectory. Skipping unwrap/center.")
        except Exception as e:
            print(f"[BioViz Notice] Skipping PBC transformations: {e}")

        # Least squares backbone alignment
        ref_universe = mda.Universe(self.topology_path, self.trajectory_path)
        align.AlignTraj(self.universe, ref_universe, select="protein and backbone",
                        in_memory=True).run()

        print("[BioViz] Backbone alignment complete.")

        if output_xtc_path:
            os.makedirs(os.path.dirname(output_xtc_path), exist_ok=True)
            with mda.Writer(output_xtc_path, len(target_selection)) as W:
                for ts in self.universe.trajectory:
                    W.write(target_selection)
            print(f"[BioViz] Exported clean trajectory to {output_xtc_path}")

        return self.universe

    def get_summary(self) -> dict:
        """Returns basic metadata regarding the loaded simulation."""
        if self.universe is None:
            return {}

        return {
            "n_atoms": len(self.universe.atoms),
            "n_frames": len(self.universe.trajectory),
            "timestep_ps": self.universe.trajectory.dt,
            "total_time_ns": (len(self.universe.trajectory) * self.universe.trajectory.dt) / 1000.0,
            "filename": self.trajectory_path
        }