import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pathlib import Path
from typing import Dict, Any, List, Optional
import imageio_ffmpeg


class TrajectoryVideoRenderer:
    """
    Compiles trajectory frames, dynamic metrics, and camera states into
    presentation-ready MP4 videos with automated chart overlays for lab meetings and investor pitches.
    """

    def __init__(self, universe: mda.Universe, metrics_data: Optional[Dict[str, Any]] = None):
        self.universe = universe
        self.metrics_data = metrics_data or {}

    def render_presentation_video(
            self,
            output_mp4_path: str,
            fps: int = 30,
            highlight_residues: Optional[List[int]] = None,
            frames_per_step: int = 10  # Repeats each frame to slow down playback smoothly at 30fps
    ) -> None:
        """
        Renders a dual-panel MP4 presentation video showing the 3D molecular structure
        alongside a synchronized 2D metric timeline.
        """
        total_steps = len(self.universe.trajectory)
        print(f"[BioViz Renderer] Initializing video compilation for {total_steps} trajectory frames (interpolated)...")

        fig, (ax_struct, ax_metric) = plt.subplots(1, 2, figsize=(12, 5))

        # Setup 2D Metric Panel (e.g., RMSD Timeline)
        time_data = self.metrics_data.get("time_ps", np.arange(total_steps))
        rmsd_data = self.metrics_data.get("rmsd", np.zeros(total_steps))

        ax_metric.plot(time_data, rmsd_data, color="#2563eb", lw=2, label="Global RMSD (Å)")
        ax_metric.set_title("Conformational Stability Timeline", fontsize=11, fontweight="bold")
        ax_metric.set_xlabel("Time (ps)")
        ax_metric.set_ylabel("RMSD (Å)")
        ax_metric.grid(True, linestyle="--", alpha=0.5)

        indicator_line = ax_metric.axvline(x=time_data[0], color="#dc2626", lw=1.5, linestyle="--",
                                           label="Current Frame")
        ax_metric.legend(loc="upper left", fontsize=9)

        # Setup 3D Structural Proxy Panel (Projection slice representation)
        protein = self.universe.select_atoms("protein")
        ca_atoms = self.universe.select_atoms("protein and name CA")

        def update_frame(interpolated_idx: int):
            # Map the smooth animation frame index back to the actual trajectory frame index
            frame_idx = interpolated_idx // frames_per_step
            if frame_idx >= total_steps:
                frame_idx = total_steps - 1

            # Advance trajectory frame
            self.universe.trajectory[frame_idx]
            current_time = float(self.universe.trajectory.time)

            # Clear structure axis and plot backbone projection
            ax_struct.clear()
            coords = ca_atoms.positions
            ax_struct.scatter(coords[:, 0], coords[:, 1], c=ca_atoms.resnums, cmap="viridis", s=15, alpha=0.8)
            ax_struct.set_title(f"BioViz Simulation State | Frame {frame_idx} ({current_time:.1f} ps)", fontsize=11,
                                fontweight="bold")
            ax_struct.set_xlabel("X Coordinate (Å)")
            ax_struct.set_ylabel("Y Coordinate (Å)")
            ax_struct.grid(True, linestyle=":", alpha=0.4)

            # Highlight specific target residues if requested
            if highlight_residues:
                highlight_atoms = self.universe.select_atoms(
                    f"protein and resnum {' '.join(map(str, highlight_residues))}")
                if len(highlight_atoms) > 0:
                    h_coords = highlight_atoms.positions
                    ax_struct.scatter(h_coords[:, 0], h_coords[:, 1], color="#dc2626", s=40, label="Target Active Site")
                    ax_struct.legend(loc="upper right", fontsize=8)

            # Update metric cursor line
            if len(time_data) > frame_idx:
                indicator_line.set_xdata([time_data[frame_idx], time_data[frame_idx]])

            return ax_struct, ax_metric

        total_animation_frames = total_steps * frames_per_step

        ani = animation.FuncAnimation(
            fig,
            update_frame,
            frames=total_animation_frames,
            interval=1000 // fps,
            blit=False
        )

        # Save output video
        output_path = Path(output_mp4_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            print(f"[BioViz Renderer] Compiling MP4 using ImageIO FFmpeg engine at {fps} fps...")
            plt.rcParams['animation.ffmpeg_path'] = imageio_ffmpeg.get_ffmpeg_exe()

            writer = animation.FFMpegWriter(fps=fps, metadata=dict(artist='BioViz'), bitrate=1800)
            ani.save(str(output_path), writer=writer, dpi=150)
            print(f"[BioViz Renderer] Presentation video successfully compiled and saved to: {output_path}")
        except Exception as e:
            print(f"[BioViz Warning] FFMpegWriter failed ({e}). Saving using pillow/gif fallback...")
            fallback_path = output_path.with_suffix(".gif")
            ani.save(str(fallback_path), writer="pillow", fps=fps)
            print(f"[BioViz Renderer] Saved animation fallback as GIF to: {fallback_path}")

        plt.close(fig)