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

using namespace cv;

int main(int argc,char *argv[])
{
	////////////////////////////////////////////////////////////////////
	// Initializing variables
	////////////////////////////////////////////////////////////////////
    double anglex = 0;
    double angley = 0;
    
    int outputw = 1920;
    int outputh = 1080;
    
    std::string tempstring;
    char anglexstr[40];
    char angleystr[40];
    
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
			infile.close();
		  }

	else std::cout << "Unable to open ini file, using defaults.";
	
	
	   
	   
} // end main
