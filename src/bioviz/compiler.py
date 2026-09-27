import cv2
import glob
import os


class VideoCompiler:
    @staticmethod
    def create_video_from_frames(image_dir: str, output_video_path: str, fps: int = 30):
        """Stitches a sequence of PNG frames into an MP4 video with strict dimension matching."""
        frame_files = sorted(glob.glob(os.path.join(image_dir, "frame_*.png")))

        if not frame_files:
            raise ValueError(f"No frames found to compile in {image_dir}")

        # Read the first frame to get dimensions
        first_frame = cv2.imread(frame_files[0])
        if first_frame is None:
            raise ValueError(f"Failed to read the first frame from {frame_files[0]}")

        orig_height, orig_width, _ = first_frame.shape

        # Ensure width and height are even numbers (required by many codecs like mp4v)
        width = orig_width if orig_width % 2 == 0 else orig_width - 1
        height = orig_height if orig_height % 2 == 0 else orig_height - 1

        # Define codec and create VideoWriter object using explicit (width, height)
        fourcc = cv2.VideoWriter_fourcc(*'mp4v')
        video_writer = cv2.VideoWriter(output_video_path, fourcc, fps, (width, height))

        if not video_writer.isOpened():
            raise RuntimeError(f"OpenCV VideoWriter failed to open for path: {output_video_path}")

        for frame_file in frame_files:
            frame = cv2.imread(frame_file)
            if frame is not None:
                # Resize frame if dimensions were adjusted to be even, or to guarantee match
                if frame.shape[1] != width or frame.shape[0] != height:
                    frame = cv2.resize(frame, (width, height))

                video_writer.write(frame)
            else:
                print(f"Warning: Skipping unreadable frame file: {frame_file}")

        video_writer.release()
        print(f"Video successfully compiled at: {output_video_path}")