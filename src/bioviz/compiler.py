import cv2
import glob
import os


class VideoCompiler:
    @staticmethod
    def create_video_from_frames(image_dir: str, output_video_path: str, fps: int = 30):
        """Stitches a sequence of PNG frames into an MP4 video."""
        frame_files = sorted(glob.glob(os.path.join(image_dir, "frame_*.png")))

        if not frame_files:
            raise ValueError(f"No frames found to compile in {image_dir}")

        # Read the first frame to get dimensions
        first_frame = cv2.imread(frame_files[0])
        height, width, _ = first_frame.shape

        # Define codec and create VideoWriter object
        fourcc = cv2.VideoWriter_fourcc(*'mp4v')
        video_writer = cv2.VideoWriter(output_video_path, fourcc, fps, (width, height))

        for frame_file in frame_files:
            frame = cv2.imread(frame_file)
            video_writer.write(frame)

        video_writer.release()
        print(f"Video successfully compiled at: {output_video_path}")