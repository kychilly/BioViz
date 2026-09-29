from bioviz import VideoCompiler


def main():
    print("Compiling existing PyMOL frames into 360-degree video...")

    # Initialize and run the video compiler
    video_status = VideoCompiler(
        frames_dir="src/bioviz/mutations/PyMOL_rendering/rotation_frames",
        output_dir="src/bioviz/mutations/PyMOL_rendering"
    ).compile_frames_to_video(output_filename="mutation_rotation_360.mp4", fps=18)

    print(video_status)


if __name__ == "__main__":
    main()