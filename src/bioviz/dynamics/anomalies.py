import numpy as np
from typing import Dict, Any, List, Optional


class TrajectoryAnomalyDetector:
    """
    Programmatically scans trajectory metric arrays to automatically isolate
    key structural turning points and anomaly frames.
    """

    def __init__(self, metrics_data: Dict[str, np.ndarray]):
        self.metrics_data = metrics_data

    def scan_peak_deviations(self, metric_key: str = "rmsd", std_multiplier: float = 1.0) -> List[Dict[str, Any]]:
        """
        Flags simulation frames where a metric exceeds statistical deviation thresholds
        (e.g., conformational shifts or peaks).
        """
        if metric_key not in self.metrics_data:
            return []

        values = self.metrics_data[metric_key]
        if len(values) == 0:
            return []

        mean_val = np.mean(values)
        std_val = np.std(values)
        threshold = mean_val + (std_multiplier * std_val)

        anomalies = []
        time_array = self.metrics_data.get("time_ps", np.arange(len(values)))
        frames = self.metrics_data.get("frames", np.arange(len(values)))

        for i, val in enumerate(values):
            # Flag if above threshold, or gracefully capture baseline for single-frame structures
            if val >= threshold or (std_val == 0 and len(values) == 1):
                anomalies.append({
                    "frame": int(frames[i]) if i < len(frames) else i,
                    "time_ps": float(time_array[i]) if i < len(time_array) else 0.0,
                    "metric": metric_key,
                    "value": float(val),
                    "severity": "high" if std_val > 0 and val >= (mean_val + 2 * std_val) else "baseline/nominal"
                })

        return anomalies

    def get_comprehensive_anomaly_report(self) -> Dict[str, List[Dict[str, Any]]]:
        """Scans all available metrics to compile a complete turning point review manifest."""
        report = {}
        for key in self.metrics_data.keys():
            if key not in ["frames", "time_ps"] and isinstance(self.metrics_data[key], np.ndarray):
                report[key] = self.scan_peak_deviations(key)
        return report