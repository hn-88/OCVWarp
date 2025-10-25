#include <iostream>
#include <opencv2/opencv.hpp>
#include <vector>
#include <string>
#include <stdio.h>

int main() {
    // --- Setup your input video ---
    cv::VideoCapture inputVideo("your_input_video.mp4");
    if (!inputVideo.isOpened()) {
        std::cerr << "Error: Could not open input video." << std::endl;
        return -1;
    }

    // --- Define output video properties ---
    int width = static_cast<int>(inputVideo.get(cv::CAP_PROP_FRAME_WIDTH));
    int height = static_cast<int>(inputVideo.get(cv::CAP_PROP_FRAME_HEIGHT));
    int fps = static_cast<int>(inputVideo.get(cv::CAP_PROP_FPS));
    const std::string outputFilename = "output.mp4";

    // --- Construct the FFmpeg command ---
    // Note: Adjust FFmpeg parameters as needed.
    // -f rawvideo: Input format is raw video frames.
    // -pixel_format bgr24: OpenCV Mat default pixel format.
    // -video_size: Dimensions of the video.
    // -framerate: Frames per second.
    // -i -: Read input from stdin.
    // -c:v hevc_nvenc: Use the NVIDIA HEVC encoder.
    // -preset slow -cq 28: Example quality settings.
    std::string ffmpeg_cmd = "ffmpeg -y -f rawvideo -pixel_format bgr24 -video_size " +
                           std::to_string(width) + "x" + std::to_string(height) +
                           " -framerate " + std::to_string(fps) +
                           " -i - -c:v hevc_nvenc -preset slow -cq 28 " +
                           outputFilename;

    // --- Open a pipe to FFmpeg ---
    // Use "w" for writing in binary mode.
    FILE* pipeout = popen(ffmpeg_cmd.c_str(), "w");
    if (!pipeout) {
        std::cerr << "Error: Could not open pipe to FFmpeg." << std::endl;
        return -1;
    }

    cv::Mat src, dst, res, dst_x, dst_y;
    // Assume res, dst_x, dst_y are initialized appropriately for your remap
    // For demonstration, we'll create empty mats of the correct size.
    res.create(height, width, CV_8UC3);
    dst_x.create(height, width, CV_32FC1);
    dst_y.create(height, width, CV_32FC1);
    // You would populate dst_x and dst_y with your warp map here.

    bool doneflag = false;
    for (;;) {
        inputVideo >> src;
        if (src.empty()) {
            doneflag = 1;
        }

        if (doneflag == 1) {
            break;
        }

        // --- Your transformations ---
        // As an example, we'll just use the source frame.
        // In your actual code, this would be your warped frame.
        // remap(res, dst, dst_x, dst_y, cv::INTER_LINEAR, cv::BORDER_CONSTANT, cv::Scalar(0, 0, 0));
        dst = src.clone(); // In your case, this will be the result of remap.


        // --- Write the frame to the FFmpeg pipe ---
        // Ensure the Mat is continuous in memory
        if (!dst.isContinuous()) {
            dst = dst.clone();
        }
        // Write the raw pixel data
        fwrite(dst.data, 1, dst.total() * dst.elemSize(), pipeout);

    } // end for(;;) loop

    // --- Close the pipe ---
    // This will signal to FFmpeg that there is no more data and it should finalize the video file.
    pclose(pipeout);

    inputVideo.release();
    std::cout << "Video processing finished." << std::endl;

    return 0;
}