#ifdef _WIN64
#include "stdafx.h"
// anything before a precompiled header is ignored, 
// so no endif here! add #endif to compile on __unix__ !
#endif
#ifdef _WIN64
#endif


/*
 * OCVWarp.cpp
 * 
 * Warps video files using the OpenCV framework. 
 * 
 * first commit:
 * Hari Nandakumar
 * 25 Jan 2020
 * 
 * 
 */

//#define _WIN64
//#define __unix__

// references 
// http://paulbourke.net/geometry/transformationprojection/
// http://www.fmwconcepts.com/imagemagick/fisheye2pano/index.php
// http://www.fmwconcepts.com/imagemagick/pano2fisheye/index.php
// 
// https://docs.opencv.org/3.4/d8/dfe/classcv_1_1VideoCapture.html
// https://docs.opencv.org/3.4/d7/d9e/tutorial_video_write.html
// https://docs.opencv.org/3.4.9/d1/da0/tutorial_remap.html
// https://stackoverflow.com/questions/60221/how-to-animate-the-command-line

#include <stdio.h>
#include <stdlib.h>

#ifdef __unix__
#include <unistd.h>
#include <libqhy/qhyccd.h>
#endif

#include <string.h>

#include <time.h>
#include <sys/stat.h>
// this is for mkdir

#include <opencv2/opencv.hpp>
#include "tinyfiledialogs.h"

using namespace cv;

void update_map( double anglex, double angley, Mat &map_x, Mat &map_y )
{
	// compute destination radius and theta 
	float rd = sqrt(map_x.cols^2+map_x.rows^2);
	/*
	# set theta so original is north rather than east
	#theta="theta=atan2(-xd,yd);"
	theta="theta=atan2(yd,xd);"


	# convert radius to phiang according to fisheye mode
	if [ "$projection" = "linear" ]; then
		# destination output diameter (dimensions) corresponds to 180 deg = pi (fov); angle is proportional to radius
		rad_per_px=`convert xc: -format "%[fx:pi/$hh]" info:`
		phiang="phiang=$rad_per_px*rd;"
	elif [ "$projection" = "stereographic" ]; then
		phiang="phiang=2*atan(2*rd/$hh);"
	elif [ "$projection" = "orthographic" ]; then
		phiang="phiang=asin(2*rd/$hh);"
	fi

	# convert theta to source (input) xs and phi to source ys
	# -rotate 90 aligns theta=0 with north and is faster than including in theta computation
	# y corresponds to h-phi, so that bottom of the input is center of output
	xs="xs=$ww2+theta*$px_per_theta;"
	ys="ys=$hh-phiang*$px_per_phi;"
	* */
	
    for( int i = 0; i < map_x.rows; i++ )
    {
        for( int j = 0; j < map_x.cols; j++ )
        {
            
                map_x.at<float>(i, j) = (float)j;
                map_y.at<float>(i, j) = (float)(map_x.rows - i);
                
                
        }
    }
    
}

int main(int argc,char *argv[])
{
	////////////////////////////////////////////////////////////////////
	// Initializing variables
	////////////////////////////////////////////////////////////////////
	bool doneflag = 0, interactivemode = 0;
    double anglex = 0;
    double angley = 0;
    
    int outputw = 1920;
    int outputh = 1080;
    
    std::string tempstring;
    char anglexstr[40];
    char angleystr[40];
    const bool askOutputType = argv[3][0] =='N';  // If false it will use the inputs codec type
    
    std::ifstream infile("OCVWarp.ini");
    
    int ind = 1;
    // inputs from ini file
    if (infile.is_open())
		  {
			
			infile >> tempstring;
			infile >> tempstring;
			infile >> tempstring;
			// first three lines of ini file are comments
			infile >> anglexstr;
			infile >> tempstring;
			infile >> angleystr;
			infile >> tempstring;
			infile >> outputw;
			infile >> tempstring;
			infile >> outputh;
			infile >> tempstring;
			infile >> interactivemode;
			infile.close();
			
			anglex = atof(anglexstr);
			angley = atof(angleystr);
		  }

	else std::cout << "Unable to open ini file, using defaults." << std::endl;
	
	namedWindow("In", 0); // 0 = WINDOW_NORMAL
	moveWindow("In", 0, 0);
	
	char const * FilterPatterns[2] =  { "*.avi","*.*" };
	char const * OpenFileName = tinyfd_openFileDialog(
		"Open a video file",
		"",
		2,
		FilterPatterns,
		NULL,
		0);

	if (! OpenFileName)
	{
		tinyfd_messageBox(
			"Error",
			"No file chosen. ",
			"ok",
			"error",
			1);
		return 1 ;
	}
	
	// reference:
	// https://docs.opencv.org/3.4/d7/d9e/tutorial_video_write.html
	
	VideoCapture inputVideo(OpenFileName);              // Open input
	if (!inputVideo.isOpened())
    {
        std::cout  << "Could not open the input video: " << OpenFileName << std::endl;
        return -1;
    }
     
    std::string OpenFileNamestr = OpenFileName;    
    std::string::size_type pAt = OpenFileNamestr.find_last_of('.');                  // Find extension point
    const std::string NAME = OpenFileNamestr.substr(0, pAt) + "F" + ".avi";   // Form the new name with container
    int ex = static_cast<int>(inputVideo.get(CAP_PROP_FOURCC));     // Get Codec Type- Int form
    // Transform from int to char via Bitwise operators
    char EXT[] = {(char)(ex & 0XFF) , (char)((ex & 0XFF00) >> 8),(char)((ex & 0XFF0000) >> 16),(char)((ex & 0XFF000000) >> 24), 0};
    Size S = Size((int) inputVideo.get(CAP_PROP_FRAME_WIDTH),    // Acquire input size
                  (int) inputVideo.get(CAP_PROP_FRAME_HEIGHT));
    Size Sout = Size(outputw,outputh);            
    VideoWriter outputVideo;                                        // Open the output
    if (askOutputType)
        outputVideo.open(NAME, ex=-1, inputVideo.get(CAP_PROP_FPS), Sout, true);
    else
        outputVideo.open(NAME, ex, inputVideo.get(CAP_PROP_FPS), Sout, true);
    if (!outputVideo.isOpened())
    {
        std::cout  << "Could not open the output video for write: " << OpenFileName << std::endl;
        return -1;
    }
    std::cout << "Input frame resolution: Width=" << S.width << "  Height=" << S.height
         << " of nr#: " << inputVideo.get(CAP_PROP_FRAME_COUNT) << std::endl;
    std::cout << "Input codec type: " << EXT << std::endl;
    int channel = 2; // Select the channel to save
    int key;
    unsigned long long framenum = 0;
    switch(argv[2][0])
    {
    case 'R' : channel = 2; break;
    case 'G' : channel = 1; break;
    case 'B' : channel = 0; break;
    }
    Mat src, res, dstr;
    std::vector<Mat> spl;
    Mat dst(Sout, CV_8UC3); // S = src.size, and src.type = CV_8UC3
    Mat map_x(Sout, CV_32FC1);
    Mat map_y(Sout, CV_32FC1);
    
    update_map(anglex, angley, map_x, map_y);
    
    for(;;) //Show the image captured in the window and repeat
    {
        inputVideo >> src;              // read
        if (src.empty()) break;         // check if at end
        //imshow("In",src);
        key = waitKey(10);
        if(interactivemode)
			update_map(anglex, angley, map_x, map_y);
        //remap( src, dst, map_x, map_y, INTER_LINEAR, BORDER_CONSTANT, Scalar(0, 0, 0) );
        //resize(src, res, Size(2048,2048), 0, 0, INTER_AREA);
        // this resize results in a 360 polar instead of 180 polar
        res=src(Range::all(), Range(1024, 3072));
         
        Point2f center( (float)res.cols / 2, (float)res.rows / 2 );
        double M = 1024;
        linearPolar(res,dst, center, M, INTER_LINEAR + WARP_FILL_OUTLIERS);
        
        std::cout << "\x1B[2K"; // Erase the entire current line.
        std::cout << "\x1B[0E"; // Move to the beginning of the current line.
        std::cout << "Frame number: " << framenum++ << " Size of res is " << res.size() << std::flush;
        rotate(dst, dstr, ROTATE_90_CLOCKWISE);
        
       //outputVideo.write(res); //save or
       outputVideo << dstr;
       
       switch (key)
				{

				case 27: //ESC key
				case 'x':
				case 'X':
					doneflag = 1;
					break;

				case 'u':
				case '+':
				case '=':	// increase angley
					angley = angley + 1.0;
					break;
					
				case 'm':
				case '-':
				case '_':	// decrease angley
					angley = angley - 1.0;
					break;
					
				case 'k':
				case '}':
				case ']':	// increase anglex
					anglex = anglex + 1.0;
					break;
					
				case 'h':
				case '{':
				case '[':	// decrease anglex
					anglex = anglex - 1.0;
					break;
				
				case 'U':
					// increase angley
					angley = angley + 10.0;
					break;
					
				case 'M':
					// decrease angley
					angley = angley - 10.0;
					break;
					
				case 'K':
					// increase anglex
					anglex = anglex + 10.0;
					break;
					
				case 'H':
					// decrease anglex
					anglex = anglex - 10.0;
					break;	
					
				default:
					break;
				
				}
				
		if (doneflag == 1)
		{
			break;
		}
    } // end for(;;) loop
    
    std::cout << std::endl << "Finished writing" << std::endl;
    return 0;
	   
	   
} // end main
