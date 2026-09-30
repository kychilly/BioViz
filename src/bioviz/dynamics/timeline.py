import MDAnalysis as mda
import numpy as np
from typing import Dict, Any, Optional


class TrajectoryTimelineController:

    def __init__(self, universe: mda.Universe, metrics_data: Optional[Dict[str, Any]] = None):
        self.universe = universe
        self.metrics_data = metrics_data or {}
        self.n_frames = len(self.universe.trajectory)

        # fast lookup maps
        self._time_array = self._extract_time_array()
        self._frame_to_time_map = {i: float(t) for i, t in enumerate(self._time_array)}
        self._time_to_frame_map = {float(t): i for i, t in enumerate(self._time_array)}

    def _extract_time_array(self) -> np.ndarray:
        if "time_ps" in self.metrics_data:
            return self.metrics_data["time_ps"]
        return np.array([ts.time for ts in self.universe.trajectory])

    def snap_to_frame(self, frame_idx: int) -> Dict[str, Any]:
        if not (0 <= frame_idx < self.n_frames):
            raise IndexError(f"Frame index {frame_idx} is out of bounds for trajectory with {self.n_frames} frames.")

        self.universe.trajectory[frame_idx]

        current_time = self._frame_to_time_map.get(frame_idx, self.universe.trajectory.time)

        frame_metrics = {}
        for key, values in self.metrics_data.items() if hasattr(self.metrics_data, 'items') else []:
            if isinstance(values, np.ndarray) and len(values) == self.n_frames:
                frame_metrics[key] = float(values[frame_idx])

        return {
            "target_frame": frame_idx,
            "time_ps": current_time,
            "metrics_at_frame": frame_metrics,
            "status": "synced"
        }

    def snap_to_time(self, target_time_ps: float) -> Dict[str, Any]:
        # Closest frame in picoseconds
        closest_frame = int(np.argmin(np.abs(self._time_array - target_time_ps)))
        return self.snap_to_frame(closest_frame)

    get_sync_payload = snap_to_frame