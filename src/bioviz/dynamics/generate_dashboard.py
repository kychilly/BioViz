import os
from pathlib import Path
import MDAnalysis as mda
from MDAnalysisTests.datafiles import GRO, XTC
import matplotlib.pyplot as plt
import numpy as np

from bioviz.dynamics.metrics import TrajectoryMetricsEngine
from bioviz.dynamics.anomalies import TrajectoryAnomalyDetector
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis


def generate_analytical_dashboard():
    print("[BioViz] Loading simulation dataset for analytical dashboard...")
    uni = mda.Universe(GRO, XTC)


    metrics_engine = TrajectoryMetricsEngine(uni)
    rmsd_res = metrics_engine.compute_rmsd()
    rmsf_res = metrics_engine.compute_rmsf()

    print("[BioViz Metrics] Analyzing hydrogen bond network dynamics...")
    hb = HydrogenBondAnalysis(
        universe=uni,
        donors_sel="protein and (name N O)",
        hydrogens_sel="protein and name H*",
        acceptors_sel="protein and (name O or name N)",
        d_a_cutoff=3.0,
        d_h_a_angle_cutoff=150.0
    )
    hb.run()

    n_frames = len(uni.trajectory)
    hbond_counts = np.zeros(n_frames)
    if hasattr(hb.results, 'hbonds') and hb.results.hbonds is not None and len(hb.results.hbonds) > 0:
        frame_indices = hb.results.hbonds[:, 0].astype(int)
        counts = np.bincount(frame_indices, minlength=n_frames)
        hbond_counts[:len(counts)] = counts

    anomaly_detector = TrajectoryAnomalyDetector(metrics_data=rmsd_res)
    anomalies = anomaly_detector.scan_peak_deviations("rmsd", std_multiplier=1.0)


    fig = plt.figure(figsize=(14, 10), constrained_layout=True)
    gs = fig.add_gridspec(nrows=2, ncols=2)

    ax_rmsd = fig.add_subplot(gs[0, :])
    ax_rmsf = fig.add_subplot(gs[1, 0])
    ax_hbond = fig.add_subplot(gs[1, 1])

    time_data = rmsd_res["time_ps"]
    rmsd_values = rmsd_res["rmsd"]
    ax_rmsd.plot(time_data, rmsd_values, color="#2563eb", lw=2, label="Global RMSD (Å)")

    if anomalies:
        anomaly_times = [a["time_ps"] for a in anomalies]
        anomaly_vals = [a["value"] for a in anomalies]
        ax_rmsd.scatter(anomaly_times, anomaly_vals, color="#dc2626", s=60, zorder=5, label="Flagged Structural Shift")

    ax_rmsd.set_title("Global Conformational Stability & Anomaly Peaks (RMSD)", fontsize=12, fontweight="bold")
    ax_rmsd.set_xlabel("Simulation Time (ps)", fontsize=10)
    ax_rmsd.set_ylabel("RMSD (Å)", fontsize=10)
    ax_rmsd.grid(True, linestyle="--", alpha=0.5)
    ax_rmsd.legend(loc="upper left", fontsize=9)

    resnums = rmsf_res["resnums"]
    rmsf_values = rmsf_res["rmsf"]
    ax_rmsf.plot(resnums, rmsf_values, color="#059669", lw=1.5)
    ax_rmsf.fill_between(resnums, rmsf_values, color="#059669", alpha=0.2)
    ax_rmsf.set_title("Local Residue Flexibility (RMSF Profile)", fontsize=12, fontweight="bold")
    ax_rmsf.set_xlabel("Residue Number", fontsize=10)
    ax_rmsf.set_ylabel("Fluctuation (Å)", fontsize=10)
    ax_rmsf.grid(True, linestyle="--", alpha=0.5)

    # --- Panel C: Hydrogen Bond Network Dynamics ---
    ax_hbond.plot(time_data, hbond_counts, color="#7c3aed", lw=1.5)
    ax_hbond.fill_between(time_data, hbond_counts, color="#7c3aed", alpha=0.2)
    ax_hbond.set_title("Hydrogen Bond Network Dynamics", fontsize=12, fontweight="bold")
    ax_hbond.set_xlabel("Simulation Time (ps)", fontsize=10)
    ax_hbond.set_ylabel("Active H-Bonds", fontsize=10)
    ax_hbond.grid(True, linestyle="--", alpha=0.5)

    output_dir = Path(__file__).parent / "sample_data"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "bioviz_analytical_dashboard.png"

    plt.savefig(str(output_path), dpi=300)
    plt.close(fig)
    print(f"\n Analytical dashboard visual successfully generated and saved at: {output_path.absolute()}")


if __name__ == "__main__":
    generate_analytical_dashboard()