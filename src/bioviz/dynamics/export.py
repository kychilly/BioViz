import MDAnalysis as mda
import numpy as np
import json
from typing import Dict, Any, Optional


class TrajectoryExportEngine:
    """
    Serializes trajectory frames, quantitative metrics, and anomaly manifests
    into lightweight JSON payloads optimized for web-based WebGL/NGLView renderers.
    """

    def __init__(self, universe: mda.Universe, metrics_data: Optional[Dict[str, Any]] = None):
        self.universe = universe
        self.metrics_data = metrics_data or {}

    def get_frame_coordinates_json(self, frame_idx: int) -> str:
        """
        Extracts lightweight atomic positions (x, y, z) for a specific frame
        and serializes them into a JSON string for real-time WebGL streaming.
        """
        if not (0 <= frame_idx < len(self.universe.trajectory)):
            raise IndexError(f"Frame index {frame_idx} out of bounds.")

        self.universe.trajectory[frame_idx]
        atoms = self.universe.select_atoms("protein")

        payload = {
            "frame": frame_idx,
            "time_ps": float(self.universe.trajectory.time),
            "n_atoms": len(atoms),
            "coordinates": atoms.positions.tolist(),  # [x, y, z] float arrays
            "resnames": atoms.resnames.tolist(),
            "resnums": atoms.resnums.tolist(),
            "names": atoms.names.tolist()
        }
        return json.dumps(payload)

    def export_full_visualization_manifest(self, output_path: str, anomalies: Optional[Dict[str, Any]] = None) -> None:
        """
        Exports a complete project state manifest containing time-series graphs,
        statistical metrics, and flagged turning points into a single frontend-ready JSON file.
        """
        # Convert numpy arrays to lists for JSON compliance
        serialized_metrics = {}
        for k, v in self.metrics_data.items():
            if isinstance(v, np.ndarray):
                serialized_metrics[k] = v.tolist()
            else:
                serialized_metrics[k] = v

        manifest = {
            "system_summary": {
                "n_atoms": len(self.universe.atoms),
                "n_frames": len(self.universe.trajectory),
                "total_time_ps": float(self.universe.trajectory.totaltime)
            },
            "metrics": serialized_metrics,
            "anomalies": anomalies or {}
        }

        with open(output_path, "w") as f:
            json.dump(manifest, f, indent=2)

        print(f"[BioViz Export] Successfully serialized full visualization manifest to: {output_path}")