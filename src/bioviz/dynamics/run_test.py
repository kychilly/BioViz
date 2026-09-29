from pathlib import Path
from bioviz.dynamics.ingestion import TrajectoryIngestionEngine
from bioviz.dynamics.metrics import TrajectoryMetricsEngine

base_dir = Path(__file__).parent

# 1. Ingest and clean trajectory using your existing pipeline
engine = TrajectoryIngestionEngine(
    topology_path=str(base_dir / "sample_data" / "4ake_apo.pdb"),
    trajectory_path=str(base_dir / "sample_data" / "4ake_apo.pdb") # or multi-frame XTC/DCD
)
uni = engine.preprocess_trajectory()

# 2. Compute quantitative metrics via backend engine
metrics_engine = TrajectoryMetricsEngine(uni)

# Calculate RMSD
rmsd_data = metrics_engine.compute_rmsd()
print(f"RMSD Calculated. Mean RMSD: {rmsd_data['rmsd'].mean():.3f} Å")

# Calculate RMSF
rmsf_data = metrics_engine.compute_rmsf()
print(f"RMSF Calculated for {len(rmsf_data['resnums'])} residues.")

# Calculate Hydrogen Bonds
hbond_data = metrics_engine.compute_hydrogen_bonds()
print(f"Hydrogen Bond Analysis complete. Average H-bonds per frame: {hbond_data['average_hbonds']:.1f}")