from pathlib import Path
from MDAnalysisTests.datafiles import GRO, XTC  # Built-in multi-frame test dataset
import MDAnalysis as mda

from bioviz.dynamics.ingestion import TrajectoryIngestionEngine
from bioviz.dynamics.metrics import TrajectoryMetricsEngine
from bioviz.dynamics.renderer import TrajectoryVideoRenderer
from bioviz.dynamics.timeline import TrajectoryTimelineController
from bioviz.dynamics.anomalies import TrajectoryAnomalyDetector
from bioviz.dynamics.export import TrajectoryExportEngine

base_dir = Path(__file__).parent

# 1. Ingestion & Preprocessing Pipeline
engine = TrajectoryIngestionEngine(
    topology_path=str(base_dir / "sample_data" / "4ake_apo.pdb"),
    trajectory_path=str(base_dir / "sample_data" / "4ake_apo.pdb")
)
uni = engine.preprocess_trajectory(output_xtc_path=str(base_dir / "sample_data" / "processed_output.xtc"))

# 2. Compute Quantitative Metrics
metrics_engine = TrajectoryMetricsEngine(uni)
rmsd_results = metrics_engine.compute_rmsd()
rmsf_results = metrics_engine.compute_rmsf()

# Combine metrics into a single dictionary dictionary container
combined_metrics = {
    **rmsd_results,
    "rmsf_values": rmsf_results["rmsf"],
    "resnums": rmsf_results["resnums"]
}

# 3. Timeline Synchronization & Anomaly Scanning
timeline = TrajectoryTimelineController(uni, metrics_data=rmsd_results)
anomaly_detector = TrajectoryAnomalyDetector(metrics_data=rmsd_results)
anomaly_report = anomaly_detector.get_comprehensive_anomaly_report()

# 4. Serialization & Frontend Export Layer
export_engine = TrajectoryExportEngine(uni, metrics_data=combined_metrics)

# Export full project analytics manifest
manifest_path = base_dir / "sample_data" / "bioviz_visualization_manifest.json"
export_engine.export_full_visualization_manifest(str(manifest_path), anomalies=anomaly_report)

# Stream single frame coordinate payload (simulating WebGL viewport stream)
frame_json_stream = export_engine.get_frame_coordinates_json(0)
print(f"[BioViz Export] Streamed Frame 0 WebGL coordinate payload ({len(frame_json_stream)} characters).")

print("NEW LINE ------------------------------------------------------")
print("[BioViz] Loading multi-frame simulation dataset (GRO/XTC)...")
uni = mda.Universe(GRO, XTC)
print(f"[BioViz] Loaded trajectory successfully! Total frames: {len(uni.trajectory)}")

# 2. Compute Quantitative Metrics across all frames
metrics_engine = TrajectoryMetricsEngine(uni)
rmsd_results = metrics_engine.compute_rmsd()
rmsf_results = metrics_engine.compute_rmsf()

combined_metrics = {
    **rmsd_results,
    "rmsf_values": rmsf_results["rmsf"],
    "resnums": rmsf_results["resnums"]
}

# 3. Timeline Synchronization & Anomaly Scanning
timeline = TrajectoryTimelineController(uni, metrics_data=rmsd_results)
anomaly_detector = TrajectoryAnomalyDetector(metrics_data=rmsd_results)
anomaly_report = anomaly_detector.get_comprehensive_anomaly_report()

# 4. Serialization & Frontend Export Layer
export_engine = TrajectoryExportEngine(uni, metrics_data=combined_metrics)
manifest_path = base_dir / "sample_data" / "bioviz_visualization_manifest.json"
export_engine.export_full_visualization_manifest(str(manifest_path), anomalies=anomaly_report)

# 5. Render Presentation Video (.mp4) with active-site highlighting
renderer = TrajectoryVideoRenderer(uni, metrics_data=rmsd_results)
video_output_path = base_dir / "sample_data" / "bioviz_presentation_reel.mp4"

renderer.render_presentation_video(
    output_mp4_path=str(video_output_path),
    fps=30,
    highlight_residues=[10, 11, 12, 13, 14, 15]  # Active-site highlights for presentation
)