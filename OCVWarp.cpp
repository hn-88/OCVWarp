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
 * Appends F to the filename and saves as default codec (DIVX avi) in the same folder.
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
// equations in figure at http://paulbourke.net/dome/dualfish2sphere/
// http://paulbourke.net/dome/dualfish2sphere/diagram.pdf

// http://www.fmwconcepts.com/imagemagick/fisheye2pano/index.php
// http://www.fmwconcepts.com/imagemagick/pano2fisheye/index.php
// 
// https://docs.opencv.org/3.4/d8/dfe/classcv_1_1VideoCapture.html
// https://docs.opencv.org/3.4/d7/d9e/tutorial_video_write.html
// https://docs.opencv.org/3.4.9/d1/da0/tutorial_remap.html
// https://stackoverflow.com/questions/60221/how-to-animate-the-command-line


// Pertinent equations from pano2fisheye:
// fov=180 for fisheye
// fov=2*phimax or phimax=fov/2
// note rmax=N/2; N=height of input
// linear: r=f*phi; f=rmax/phimax; f=(N/2)/((fov/2)*(pi/180))=N*180/(fov*pi)
// substitute fov=180
// linear: f=N/pi
// linear: phi=r*pi/N

// https://stackoverflow.com/questions/46883320/conversion-from-dual-fisheye-coordinates-to-equirectangular-coordinates
// taking Paul's page as ref, http://paulbourke.net/dome/dualfish2sphere/diagram.pdf
/* // 2D fisheye to 3D vector
phi = r * aperture / 2
theta = atan2(y, x)

// 3D vector to longitude/latitude
longitude = atan2(Py, Px)
latitude = atan2(Pz, (Px^2 + Py^2)^(0.5))

// 3D vector to 2D equirectangular
x = longitude / PI
y = 2 * latitude / PI
* ***/
/*
 * https://groups.google.com/forum/#!topic/hugin-ptx/wB-4LJHH5QI
 * panotools code
 * */

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

void update_map( double anglex, double angley, Mat &map_x, Mat &map_y, int transformtype )
{
	// explanation comments are most verbose in the last 
	// default (transformtype == 0) section
	
	if (transformtype == 1)	// Equirectangular 360 to 180 degree fisheye
	{
		int xcd = floor(map_x.cols/2) - 1 + anglex;
		int ycd = floor(map_x.rows/2) - 1 + angley;
		int xd, yd;
		float px_per_theta = map_x.cols * 2 / (2*CV_PI); 	// src width = map_x.cols * 2
		float px_per_phi   = map_x.rows / CV_PI;			// src height = PI for equirect 360
		float rad_per_px = CV_PI / map_x.rows;
		float rd, theta, phiang;
		for ( int i = 0; i < map_x.rows; i++ ) // here, i is for y and j is for x
			{
				for ( int j = 0; j < map_x.cols; j++ )
				{
					xd = j - xcd;
					yd = i - ycd;
					if (xd == 0 && yd == 0)
					{
						theta = 0;
						rd = 0;
					}
					else
					{
						theta = atan2(xd,yd); // this sets orig to north
						rd = sqrt(float(xd*xd + yd*yd));
					}
					
					phiang = rad_per_px * rd;
					
					map_x.at<float>(i, j) = (float)round((map_x.cols) + theta * px_per_theta);
					// halfway point of src = map_x.cols
					map_y.at<float>(i, j) = phiang * px_per_phi;
					
				 } // for j
				   
			} // for i
			
	}
	else
	//if (transformtype == 0) is the default // Equirectangular 360 to 360 degree fisheye
	{
		
			// set destination (output) centers
			int xcd = floor(map_x.cols/2) - 1;
			int ycd = floor(map_x.rows/2) - 1;
			int xd, yd;
			//define destination (output) coordinates center relative xd,yd
			// "xd= x - xcd;"
			// "yd= y - ycd;"

			// compute input pixels per angle in radians
			// theta ranges from -180 to 180 = 360 = 2*pi
			// phi ranges from 0 to 90 = pi/2
			float px_per_theta = map_x.cols / (2*CV_PI);
			float px_per_phi   = map_x.rows / (CV_PI/2);
			// compute destination radius and theta 
			float rd; // = sqrt(x^2+y^2);
			
			// set theta so original is north rather than east
			float theta; //= atan2(y,x);
			
			// convert radius to phiang according to fisheye mode
			//if projection is linear then
			//	 destination output diameter (dimensions) corresponds to 180 deg = pi (fov); angle is proportional to radius
			float rad_per_px = CV_PI / map_x.rows;
			float phiang;     // = rad_per_px * rd;
			

			// convert theta to source (input) xs and phi to source ys
			// -rotate 90 aligns theta=0 with north and is faster than including in theta computation
			// y corresponds to h-phi, so that bottom of the input is center of output
			// xs = width + theta * px_per_theta;
			// ys = height - phiang * px_per_phi;
			
			
			for ( int i = 0; i < map_x.rows; i++ ) // here, i is for y and j is for x
			{
				for ( int j = 0; j < map_x.cols; j++ )
				{
					xd = j - xcd;
					yd = i - ycd;
					if (xd == 0 && yd == 0)
					{
						theta = 0;
						rd = 0;
					}
					else
					{
						//theta = atan2(float(yd),float(xd)); // this sets orig to east
						// so America, at left of globe, becomes centred
						theta = atan2(xd,yd); // this sets orig to north
						// makes the fisheye left/right flipped if atan2(-xd,yd)
						// so that Africa is centred.
						rd = sqrt(float(xd*xd + yd*yd));
					}
					
					phiang = rad_per_px * rd;
					
					map_x.at<float>(i, j) = (float)round((map_x.cols/2) + theta * px_per_theta);
					
					//map_y.at<float>(i, j) = (float)round((map_x.rows) - phiang * px_per_phi);
					// this above makes the south pole the centre.
					
					map_y.at<float>(i, j) = phiang * px_per_phi;
					// this above makes the north pole the centre of the fisheye
					
					 
				   // the following test mapping just makes the src upside down in dst
				   // map_x.at<float>(i, j) = (float)j;
				   // map_y.at<float>(i, j) = (float)( i); 
				   
				 } // for j
				   
			} // for i
				
            
     
     } // end of if transformtype == 0
    
    
// debug
    /*
    std::cout << "map_x -> " << std::endl;
    
    for ( int i = 0; i < map_x.rows; i+=100 ) // here, i is for y and j is for x
    {
        for ( int j = 0; j < map_x.cols; j+=100 )
        {
			std::cout << map_x.at<float>(i, j) << " " ;
		}
		std::cout << std::endl;
	}
	
	std::cout << "map_y -> " << std::endl;
    
    for ( int i = 0; i < map_x.rows; i+=100 ) // here, i is for y and j is for x
    {
        for ( int j = 0; j < map_x.cols; j+=100 )
        {
			std::cout << map_y.at<float>(i, j) << " " ;
		}
		std::cout << std::endl;
	}
	* */    
    
} // end function updatemap

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
    
    int transformtype = 0;
    // 0 = Equirectangular to 360 degree fisheye
    // 1 = Equirectangular to 180 degree fisheye
    
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
			infile >> tempstring;
			infile >> transformtype;
			infile.close();
			
			anglex = atof(anglexstr);
			angley = atof(angleystr);
		  }

	else std::cout << "Unable to open ini file, using defaults." << std::endl;
	
	namedWindow("Display", WINDOW_NORMAL | WINDOW_KEEPRATIO); // 0 = WINDOW_NORMAL
	resizeWindow("Display", 640, 640); 
	moveWindow("Display", 0, 0);
	
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
     
    int  fps, key;
	int t_start, t_end;
    unsigned long long framenum = 0;
     
    Mat src, res;
    std::vector<Mat> spl;
    Mat dst(Sout, CV_8UC3); // S = src.size, and src.type = CV_8UC3
    Mat map_x(Sout, CV_32FC1);
    Mat map_y(Sout, CV_32FC1);
    
    update_map(anglex, angley, map_x, map_y, transformtype);
    t_start = time(NULL);
	fps = 0;
	
    for(;;) //Show the image captured in the window and repeat
    {
        inputVideo >> src;              // read
        if (src.empty()) break;         // check if at end
        //imshow("Display",src);
        key = waitKey(10);
        
        if(interactivemode)
			update_map(anglex, angley, map_x, map_y, transformtype);
		
		switch (transformtype)
				{

				case 1: // Equirect to 180 fisheye
					resize( src, res, Size(outputw*2, outputh), 0, 0, INTER_AREA);
					break;
				
				default:	
				case 0: // Equirect to 360 fisheye
					resize( src, res, Size(outputw, outputh), 0, 0, INTER_AREA);
					break;
					
				}
		
        remap( res, dst, map_x, map_y, INTER_LINEAR, BORDER_CONSTANT, Scalar(0, 0, 0) );
        
        imshow("Display", dst);
        //std::cout << "\x1B[2K"; // Erase the entire current line.
        std::cout << "\x1B[0E"; // Move to the beginning of the current line.
        fps++;
        t_end = time(NULL);
		if (t_end - t_start >= 5)
		{
			std::cout << "Frame number: " << framenum++ << "x:" << anglex << "y:" << angley << " fps: " << fps/5 << std::flush;
			t_start = time(NULL);
			fps = 0;
		}
		//else
        std::cout << "Frame number: " << framenum++ << "x:" << anglex << "y:" << angley << std::flush;
        
        
       //outputVideo.write(res); //save or
       outputVideo << dst;
       
       
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
