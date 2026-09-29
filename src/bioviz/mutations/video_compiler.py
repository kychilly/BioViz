import os
import cv2


class VideoCompiler:
    """Takes existing pre-rendered PNG frame sequences and compiles them into an MP4 video using OpenCV."""

    def __init__(self, frames_dir: str = "src/bioviz/mutations/PyMOL_rendering/rotation_frames", output_dir: str = "src/bioviz/mutations/PyMOL_rendering"):
        self.frames_dir = os.path.abspath(frames_dir)
        self.output_dir = os.path.abspath(output_dir)

        # Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)

    def compile_frames_to_video(self, output_filename: str = "mutation_rotation_360.mp4", fps: int = 24) -> str:
        """
        Scans the frames directory for sequential PNG files, reads their dimensions,
        and writes them out into an MP4 video file.
        """
        if not os.path.exists(self.frames_dir):
            return f"Error: Frames directory not found at {self.frames_dir}"

        # Get all PNG images and sort them alphabetically so they align with frame naming (frame_000, frame_001, etc.)
        images = sorted([img for img in os.listdir(self.frames_dir) if img.endswith(".png")])
        if not images:
            return f"Error: No PNG frames found inside {self.frames_dir}"

        output_video_path = os.path.join(self.output_dir, output_filename)
        print(f"Found {len(images)} frames. Compiling into video at {fps} FPS...")

        # Read the first frame to determine video frame width and height dynamically
        first_frame_path = os.path.join(self.frames_dir, images[0])
        sample_img = cv2.imread(first_frame_path)
        if sample_img is None:
            return f"Error: Could not read sample frame at {first_frame_path}"

        height, width, _ = sample_img.shape

        # Initialize OpenCV video writer (using mp4v codec)
        fourcc = cv2.VideoWriter_fourcc(*'mp4v')
        video_writer = cv2.VideoWriter(output_video_path, fourcc, fps, (width, height))

        for image in images:
            img_path = os.path.join(self.frames_dir, image)
            frame = cv2.imread(img_path)
            if frame is not None:
                video_writer.write(frame)

        video_writer.release()
        return f"Video successfully compiled and saved at: {output_video_path}"