#ifdef _WIN64
#include "windows.h"
#endif



/*
 * modified from OCVWarp.cpp
 * 
 * Outputs the list of fourccs from opencv.
 * 
 * first commit:
 * Hari Nandakumar
 * 2 Apr 2020
 * 
 * 
 */

//#define _WIN64
//#define __unix__

// references 
// https://docs.opencv.org/3.4/d8/dfe/classcv_1_1VideoCapture.html
// https://docs.opencv.org/3.4/d7/d9e/tutorial_video_write.html
// https://docs.opencv.org/3.4.9/d1/da0/tutorial_remap.html
// https://stackoverflow.com/questions/60221/how-to-animate-the-command-line

#include <stdio.h>
#include <stdlib.h>

#ifdef __unix__
#include <unistd.h>
#endif

#include <string.h>

#include <time.h>
//#include <sys/stat.h>
// this is for mkdir

#include <opencv2/opencv.hpp>
#include "tinyfiledialogs.h"

#define CV_PI   3.1415926535897932384626433832795

using namespace cv;


int main(int argc,char *argv[])
{
	
	// reference:
	// https://docs.opencv.org/3.4/d7/d9e/tutorial_video_write.html
	bool askOutputType = 1;
	std::string OpenFileNamestr = "input.mp4";    
    std::string::size_type pAt = OpenFileNamestr.find_last_of('.');                  // Find extension point
    const std::string NAME = OpenFileNamestr.substr(0, pAt) + "F" + ".avi";   // Form the new name with container
    int ex  ; // = static_cast<int>(inputVideo.get(CAP_PROP_FOURCC));     // Get Codec Type- Int form
    // Transform from int to char via Bitwise operators
    char EXT[] = {(char)(ex & 0XFF) , (char)((ex & 0XFF00) >> 8),(char)((ex & 0XFF0000) >> 16),(char)((ex & 0XFF000000) >> 24), 0};
    Size S = Size((int) 640,    // Acquire input size
                  (int) 480);
    Size Sout = Size(640,480);            
    VideoWriter outputVideo;                                        // Open the output
    if (askOutputType)
        outputVideo.open(NAME, -1 , 25, Sout, true);
    else
        outputVideo.open(NAME, ex, 25, Sout, true);
    if (!outputVideo.isOpened())
    {
        std::cout  << "Could not open the output video for write: " << NAME << std::endl;
        return -1;
    }
       
	   
} // end main
