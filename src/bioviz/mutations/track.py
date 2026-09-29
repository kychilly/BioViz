class GenomeTrackRenderer:
    """Processes mutation density, clustering hotspots, and generates linear timeline tracks."""

    def __init__(self, sequence_length: int, window_size: int = 10):
        self.sequence_length = sequence_length
        self.window_size = window_size

    def compute_density_and_hotspots(self, annotated_variants: list) -> dict:
        """Bins variants into genomic windows to calculate density and flag clustering hotspots."""
        num_windows = (self.sequence_length // self.window_size) + 1
        windows = [0] * num_windows
        hotspot_threshold = 2  # Threshold to classify a region as a mutation hotspot

        # Map variants to their respective sliding-window bins
        for v in annotated_variants:
            pos = v["position"] - 1
            window_idx = pos // self.window_size
            if window_idx < num_windows:
                windows[window_idx] += 1

        track_data = []
        hotspots = []

        for i, count in enumerate(windows):
            start_bp = i * self.window_size + 1
            end_bp = min((i + 1) * self.window_size, self.sequence_length)

            is_hotspot = count >= hotspot_threshold
            if is_hotspot:
                hotspots.append({
                    "window": i,
                    "start_bp": start_bp,
                    "end_bp": end_bp,
                    "mutation_count": count
                })

            track_data.append({
                "window_index": i,
                "region": f"{start_bp}-{end_bp}",
                "mutation_count": count,
                "hotspot": is_hotspot
            })

        return {
            "window_size": self.window_size,
            "total_length": self.sequence_length,
            "track_bins": track_data,
            "identified_hotspots": hotspots
        }

    @staticmethod
    def render_ascii_timeline(track_result: dict) -> str:
        """Generates a scannable ASCII text timeline map for local terminal output."""
        ascii_bar = "\n--- Genomic Coordinate Track & Timeline ---\n"
        ascii_bar += "[ 5' "

        for bin_info in track_result["track_bins"]:
            count = bin_info["mutation_count"]
            if bin_info["hotspot"]:
                ascii_bar += f" 🔥[{count}] "
            elif count > 0:
                ascii_bar += f" •[{count}] "
            else:
                ascii_bar += " --- "

        ascii_bar += " 3' ]\n"
        return ascii_bar