#ifdef _WIN64
#include "windows.h"
#endif

/*
 * OCVWarp.cpp
 * 
 * Warps video files using the OpenCV framework. 
 * Documentation is at https://github.com/hn-88/OCVWarp/wiki
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
// https://stackoverflow.com/questions/11498169/dealing-with-angle-wrap-in-c-code
// https://blog.kowalczyk.info/article/j/guide-to-predefined-macros-in-c-compilers-gcc-clang-msvc-etc..html

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
#include <fstream>
#include <time.h>
//#include <sys/stat.h>
// this is for mkdir

#include <opencv2/opencv.hpp>
#include "tinyfiledialogs.h"
#define CVUI_IMPLEMENTATION
#include "cvui.h"
#define WINDOW_NAME "OCVWARP - HIT <esc> TO CLOSE"

#define CV_PI   3.1415926535897932384626433832795

using namespace cv;

// some global variables

std::string strpathtowarpfile;
Mat meshu, meshv, meshx, meshy, meshi, I;
Mat map2x, map2y;
float maxx=0, minx=0;	
float maxu=0, minu=0, maxv=0, minv=0;
// meshx is in the range [-aspectratio, aspectratio]
// we assume meshy is in the range [-1,1]
// meshu and meshv in [0,1]
std::string tempstring;

int outputw = 1920;
int outputh = 1080;

int transformtype = 0;
 	// 0 = Equirectangular to 360 degree fisheye
    	// 1 = Equirectangular to 180 degree fisheye
    
    char anglexstr[40];
    char angleystr[40];
    char anglexincrstr[40];
    char angleyincrstr[40];
    char outputfourccstr[40];	// leaving extra chars for not overflowing too easily
    char outputfpsstr[40];    

bool ReadMesh(std::string strpathtowarpfile)
{
	
	//from https://github.com/hn-88/GL_warp2Avi/blob/master/GL2AviView.cpp
	// and http://paulbourke.net/dataformats/meshwarp/
	FILE *input = NULL;

   input = fopen(strpathtowarpfile.c_str(), "r");

   /* Set rows and columns to 2 initially, as this is the size of the default mesh. */
    int dummy, rows = 2, cols = 2;

    if (input != NULL)  {
		fscanf(input, " %d %d %d ", &dummy, &cols, &rows) ;
		float x, y, u, v, l;
		//meshrows=rows;
		//meshcolumns=cols;

		meshx = Mat(Size(cols,rows), CV_32FC1);	
		meshy = Mat(Size(cols,rows), CV_32FC1);
		meshu = Mat(Size(cols,rows), CV_32FC1);
		meshv = Mat(Size(cols,rows), CV_32FC1);
		meshi = Mat(Size(cols,rows), CV_32FC1);

	     for (int r = 0; r < rows ; r++) {
             for (int c = 0; c < cols ; c++) {
		                
                fscanf(input, "%f %f %f %f %f", &x, &y, &u, &v, &l) ;   
                
                if (x<minx)
					minx = x;
				else
				if (x>maxx)
					maxx = x;
					
				if (u<minu)
					minu = u;
				else
				if (u>maxu)
					maxu = u;
					
				if (v<minv)
					minv = v;
				else
				if (v>maxv)
					maxv = v;
                
                //~ mesh[cols*r+c].x = x;
                //~ mesh[cols*r+c].y = y;
                //~ mesh[cols*r+c].u = u;
                //~ mesh[cols*r+c].v = v;
                //~ mesh[cols*r+c].i = l;
                meshx.at<float>(r,c) = x;
                meshy.at<float>(r,c) = y;
                meshu.at<float>(r,c) = u;
                meshv.at<float>(r,c) = v;
                meshi.at<float>(r,c) = l;
 
			}
			 
		}
		
	}
	else // unable to read mesh
	{
		std::cout << "Unable to read mesh data file (similar to EP_xyuv_1920.map), exiting!" << std::endl;
		exit(0);
	}
	
	return 1;
}

void update_map( double anglex, double angley, Mat &map_x, Mat &map_y, int transformtype )
{
	// explanation comments are most verbose in the last 
	// default (transformtype == 0) section
	switch (transformtype)
	{
	case 5: // 360 to 180 fisheye and then to warped
	//if (transformtype == 5)	
	{
		// create temp maps to the texture and then map from texture to output
		// this will need 2 remaps at the output side
		// and two sets of map files
		
		// the map file for the first remap has to change with change of anglex angley
		// the second one, for fisheye to warped, doesn't need to be recalculated.
		
		// so, update_map is called first to initialize map_x and map_y
		// using transformtype = 4
		
		// and later, and at all other times, only the map to the texture is updated.
		
		update_map( anglex, angley, map2x, map2y, 1 );
		
	}
		break;
		
	case 4:
	//if (transformtype == 4)	//  fisheye to warped
	{
		// similar to TGAWarp at http://paulbourke.net/dome/tgawarp/
		//
		
		Mat U, V, X, Y, IC1;
		Mat indexu, indexv, indexx, indexy, temp;
		ReadMesh(strpathtowarpfile);
		//resize(meshx, X, map_x.size(), INTER_LINEAR);
		//resize(meshy, Y, map_x.size(), INTER_LINEAR);
		//debug - changed INTER_LINEAR to INTER_LANCZOS4 and later INTER_CUBIC
		// not much of a penalty, so we leave it in.
		
		// discard the top/bottom line of U and V, since they cause
		// the bottom of the image to be repeats of the same
		//~ Mat meshub, meshvb;
		//~ meshu(cv::Rect(0,0,meshu.cols,(meshu.rows-1))).copyTo(meshub);
		//~ meshv(cv::Rect(0,0,meshv.cols,(meshv.rows-1))).copyTo(meshvb);
		// this doesn't work
		
		// for per pixel equivalence with GL_warp, the following seems to be needed
		int extrarows = map_x.rows / meshu.rows ;
		int extracols = map_x.cols / meshu.cols;
		resize(meshu, U, Size(map_x.cols+extracols, (map_x.rows+extrarows)), INTER_CUBIC);
		resize(meshv, V, Size(map_x.cols+extracols, (map_x.rows+extrarows)), INTER_CUBIC);
		
		//resize(meshu, U, map_x.size(), INTER_CUBIC);
		//resize(meshv, V, map_x.size(), INTER_CUBIC);
		resize(meshi, IC1, map_x.size(), INTER_CUBIC);
		
		// I.convertTo(I, CV_32FC3); //this doesn't work 	
		//convert to 3 channel Mat, for later multiplication
		// https://stackoverflow.com/questions/23303305/opencv-element-wise-matrix-multiplication/23310109
		Mat t[] = {IC1, IC1, IC1};
		merge(t, 3, I);
		
		// map the values which are [minx,maxx] to [0,map_x.cols-1]
		//~ temp = (map_x.cols-1)*(X - minx)/(maxx-minx);
		//~ temp.convertTo(indexx, CV_32S);		// this does the rounding to int
		
		//~ temp = (map_x.rows-1)*(Y+1)/2;	// assuming miny=-1, maxy=1
		//~ temp.convertTo(indexy, CV_32S);
		
		temp = (map_x.cols-1)*U;	// assuming minu=0, maxu=1
		temp.convertTo(indexu, CV_32S);		// this does the rounding to int
		
		temp = (map_x.rows-1)*V;	// assuming minv=0, maxv=1
		//~ temp = (map_x.rows)*(V-minv)/(maxv-minv);
		temp.convertTo(indexv, CV_32S);
		
		int linestodiscard = map_x.rows / meshu.rows / 2;
		int colstodiscard = map_x.cols / meshu.cols / 2;
		
		for ( int i = 0; i < (map_x.rows); i++ ) // here, i is for y and j is for x
			{
				for ( int j = 0; j < map_x.cols; j++ )
				{
					//~ map_x.at<float>(i, j) = (float)(j); // this just maps input to output
				    //~ map_y.at<float>(i, j) = (float)(i); 
				    
				    // in the following, we assume indexx.at<int>(i,j) = j
				    // and indexy.at<int>(i,j) = i
				    // otherwise, a mesh effect due to discontinuities in indexx and indexy.
				    	
					map_x.at<float>(i,j) = (float) indexu.at<int>(i+linestodiscard,j+colstodiscard);
					map_y.at<float>(i,j) = (float) indexv.at<int>(i+linestodiscard,j+colstodiscard);
					
				} //end for j
			} //end for i
		return;
	}
	break;
		
	case 3:
	//if (transformtype == 3)	//  fisheye to Equirectangular - dual output - using parallel projection
	{
		// int xcd = floor(map_x.cols/2) - 1 + anglex;	// this just 'pans' the view
		// int ycd = floor(map_x.rows/2) - 1 + angley;
		int xcd = floor(map_x.cols/2) - 1 ;
		int ycd = floor(map_x.rows/2) - 1 ;
		 
		float px_per_theta = map_x.cols * 2 / (2*CV_PI); 	// src width = map_x.cols * 2
		float px_per_phi   = map_x.rows / CV_PI;			// src height = PI for equirect 360
		float rad_per_px = CV_PI / map_x.rows;
		float theta;
		float longi, lat, Px, Py, Pz, R;						// X and Y are map_x and map_y
		float PxR, PyR, PzR;
		float aperture = CV_PI;	// this is the only change between type 2 & 3
		float angleyrad = -angley*CV_PI/180;	// made these minus for more intuitive feel
		float anglexrad = -anglex*CV_PI/180;
		
		
		for ( int i = 0; i < map_x.rows; i++ ) // here, i is for y and j is for x
			{
				for ( int j = 0; j < map_x.cols; j++ )
				{
					longi 	= (CV_PI    ) * (j - xcd) / (map_x.cols/2);		// longi = x.pi for 360 image
					lat	 	= (CV_PI / 2) * (i - ycd) / (map_x.rows/2);		// lat = y.pi/2
					
					Px = cos(lat)*cos(longi);
					Py = cos(lat)*sin(longi);
					Pz = sin(lat);
					
					if(angley!=0 || anglex!=0)
					{
						// cos(angleyrad), 0, sin(angleyrad), 0, 1, 0, -sin(angleyrad), 0, cos(angleyrad));
						
						PxR = Px;
						PyR = cos(angleyrad) * Py - sin(angleyrad) * Pz;
						PzR = sin(angleyrad) * Py + cos(angleyrad) * Pz;
						
						Px = cos(anglexrad) * PxR - sin(anglexrad) * PyR;
						Py = sin(anglexrad) * PxR + cos(anglexrad) * PyR;
						Pz = PzR;
					}
					
					
					if (Px == 0 && Py == 0 && Pz == 0)
						R = 0;
					else 
						R = 2 * atan2(sqrt(Px*Px + Pz*Pz), Py) / aperture; 	
					
					if (Px == 0 && Pz ==0)
						theta = 0;
					else
						theta = atan2(Pz, Px);
						
					
					// map_x.at<float>(i, j) = R * cos(theta); this maps to [-1, 1]
					//map_x.at<float>(i, j) = R * cos(theta) * map_x.cols / 2 + xcd;
					map_x.at<float>(i, j) = - Px * map_x.cols / 2 + xcd;
					
					// this gives two copies in final output, top one reasonably correct
					
					// map_y.at<float>(i, j) = R * sin(theta); this maps to [-1, 1]
					//map_y.at<float>(i, j) = R * sin(theta) * map_x.rows / 2 + ycd;
					map_y.at<float>(i, j) = Py * map_x.rows / 2 + ycd;
					
				 } // for j
				   
			} // for i
			
	}
	break;
	
	case 2:
	//if (transformtype == 2)	// 360 degree fisheye to Equirectangular 360 
	{
		// int xcd = floor(map_x.cols/2) - 1 + anglex;	// this just 'pans' the view
		// int ycd = floor(map_x.rows/2) - 1 + angley;
		int xcd = floor(map_x.cols/2) - 1 ;
		int ycd = floor(map_x.rows/2) - 1 ;
		 
		float px_per_theta = map_x.cols  / (2*CV_PI); 	//  width = map_x.cols 
		float px_per_phi   = map_x.rows / CV_PI;		//  height = PI for equirect 360
		float rad_per_px = CV_PI / map_x.rows;
		float theta;
		float longi, lat, Px, Py, Pz, R;						// X and Y are map_x and map_y
		float PxR, PyR, PzR;
		float aperture = 2*CV_PI;
		float angleyrad = -angley*CV_PI/180;	// made these minus for more intuitive feel
		float anglexrad = -anglex*CV_PI/180;
		
		for ( int i = 0; i < map_x.rows; i++ ) // here, i is for y and j is for x
			{
				for ( int j = 0; j < map_x.cols; j++ )
				{
					longi 	= (CV_PI    ) * (j - xcd) / (map_x.cols/2);		// longi = x.pi for 360 image
					lat	 	= (CV_PI / 2) * (i - ycd) / (map_x.rows/2);		// lat = y.pi/2
					
					Px = cos(lat)*cos(longi);
					Py = cos(lat)*sin(longi);
					Pz = sin(lat);
					
					if(angley!=0 || anglex!=0)
					{
						// cos(angleyrad), 0, sin(angleyrad), 0, 1, 0, -sin(angleyrad), 0, cos(angleyrad));
						
						PxR = Px;
						PyR = cos(angleyrad) * Py - sin(angleyrad) * Pz;
						PzR = sin(angleyrad) * Py + cos(angleyrad) * Pz;
						
						Px = cos(anglexrad) * PxR - sin(anglexrad) * PyR;
						Py = sin(anglexrad) * PxR + cos(anglexrad) * PyR;
						Pz = PzR;
					}
					
					
					if (Px == 0 && Py == 0 && Pz == 0)
						R = 0;
					else 
						R = 2 * atan2(sqrt(Px*Px + Py*Py), Pz) / aperture; 
						// exchanged Py and Pz from Paul's co-ords, 	
						// from Perspective projection the wrong imaging model 10.1.1.52.8827.pdf
						// http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.52.8827&rep=rep1&type=pdf
						// Or else, Africa ends up sideways, and with the far east and west streched out on top and bottom
					
					if (Px == 0 && Pz ==0)
						theta = 0;
					else
						theta = atan2(Py, Px);
						
					
					// map_x.at<float>(i, j) = R * cos(theta); this maps to [-1, 1]
					map_x.at<float>(i, j) =  R * cos(theta) * map_x.cols / 2 + xcd;
					
					// currently upside down 
					
					// map_y.at<float>(i, j) = R * sin(theta); this maps to [-1, 1]
					map_y.at<float>(i, j) =  R * sin(theta) * map_x.rows / 2 + ycd;
					
					
				 } // for j
				   
			} // for i
			
	}
	break;
	
	case 1:
	//if (transformtype == 1)	// Equirectangular 360 to 180 degree fisheye
	{
		// using the transformations at
		// http://paulbourke.net/dome/dualfish2sphere/diagram.pdf
		int xcd = floor(map_x.cols/2) - 1 ;
		int ycd = floor(map_x.rows/2) - 1 ;
		float halfcols = map_x.cols/2;
		float halfrows = map_x.rows/2;
		
		
		float longi, lat, Px, Py, Pz, theta;						// X and Y are map_x and map_y
		float xfish, yfish, rfish, phi, xequi, yequi;
		float PxR, PyR, PzR;
		float aperture = CV_PI;
		float angleyrad = -angley*CV_PI/180;	// made these minus for more intuitive feel
		float anglexrad = -anglex*CV_PI/180;
		
		//Mat inputmatrix, rotationmatrix, outputmatrix;
		// https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations
		//rotationmatrix = (Mat_<float>(3,3) << cos(angleyrad), 0, sin(angleyrad), 0, 1, 0, -sin(angleyrad), 0, cos(angleyrad)); //y
		//rotationmatrix = (Mat_<float>(3,3) << 1, 0, 0, 0, cos(angleyrad), -sin(angleyrad), 0, sin(angleyrad), cos(angleyrad)); //x
		//rotationmatrix = (Mat_<float>(3,3) << cos(angleyrad), -sin(angleyrad), 0, sin(angleyrad), cos(angleyrad), 0, 0, 0, 1); //z
		
		for ( int i = 0; i < map_x.rows; i++ ) // here, i is for y and j is for x
			{
				for ( int j = 0; j < map_x.cols; j++ )
				{
					// normalizing to [-1, 1]
					xfish = (j - xcd) / halfcols;
					yfish = (i - ycd) / halfrows;
					rfish = sqrt(xfish*xfish + yfish*yfish);
					theta = atan2(yfish, xfish);
					phi = rfish*aperture/2;
					
					// Paul's co-ords - this is suitable when phi=0 is Pz=0
					
					//Px = cos(phi)*cos(theta);
					//Py = cos(phi)*sin(theta);
					//Pz = sin(phi);
					
					// standard co-ords - this is suitable when phi=pi/2 is Pz=0
					Px = sin(phi)*cos(theta);
					Py = sin(phi)*sin(theta);
					Pz = cos(phi);
					
					if(angley!=0 || anglex!=0)
					{
						// cos(angleyrad), 0, sin(angleyrad), 0, 1, 0, -sin(angleyrad), 0, cos(angleyrad));
						
						PxR = Px;
						PyR = cos(angleyrad) * Py - sin(angleyrad) * Pz;
						PzR = sin(angleyrad) * Py + cos(angleyrad) * Pz;
						
						Px = cos(anglexrad) * PxR - sin(anglexrad) * PyR;
						Py = sin(anglexrad) * PxR + cos(anglexrad) * PyR;
						Pz = PzR;
					}
					
					
					longi 	= atan2(Py, Px);
					lat	 	= atan2(Pz,sqrt(Px*Px + Py*Py));	
					// this gives south pole centred, ie yequi goes from [-1, 0]
					// Made into north pole centred by - (minus) in the final map_y assignment
					
					xequi = longi / CV_PI;
					// this maps to [-1, 1]
					yequi = 2*lat / CV_PI;
					// this maps to [-1, 0] for south pole
					
					//if (rfish <= 1.0)		// outside that circle, let it be black
					// removed the black circle to help transformtype=5
					// avoid bottom pixels black
					{
						map_x.at<float>(i, j) =  abs(xequi * map_x.cols / 2 + xcd);
						//map_y.at<float>(i, j) =  yequi * map_x.rows / 2 + ycd;
						// this gets south pole centred view
						
						// the abs is to correct for -0.5 xequi value at longi=0
						
						map_y.at<float>(i, j) =  yequi * map_x.rows / 2 + ycd;
						//debug
						//~ if (rfish <= 1.0/500)
						//if ((longi==0)||(longi==CV_PI)||(longi==-CV_PI))
						//if (lat==0)	// since these are floats, probably doesn't work
						//~ {
							//~ std::cout << "i,j,mapx,mapy=";
							//~ std::cout << i << ", ";
							//~ std::cout << j << ", ";
							//~ std::cout << map_x.at<float>(i, j) << ", ";
							//~ std::cout << map_y.at<float>(i, j) << std::endl;
						//~ }
					}
					
				 } // for j
				   
			} // for i
			
	}
	break;
	
	case 0:
	default:
	//else
	//if (transformtype == 0) // the default // Equirectangular 360 to 360 degree fisheye
	{
			
		//////////////////////////////////////
		// the following code is similar to transformtype=1 code
		// with only the "aperture" changed to 2pi	
		int xcd = floor(map_x.cols/2) - 1 ;
		int ycd = floor(map_x.rows/2) - 1 ;
		float halfcols = map_x.cols/2;
		float halfrows = map_x.rows/2;
		
		
		float longi, lat, Px, Py, Pz, theta;						// X and Y are map_x and map_y
		float xfish, yfish, rfish, phi, xequi, yequi;
		float PxR, PyR, PzR;
		float aperture = 2*CV_PI;
		float angleyrad = -angley*CV_PI/180;	// made these minus for more intuitive feel
		float anglexrad = -anglex*CV_PI/180;
		
		//Mat inputmatrix, rotationmatrix, outputmatrix;
		// https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations
		//rotationmatrix = (Mat_<float>(3,3) << cos(angleyrad), 0, sin(angleyrad), 0, 1, 0, -sin(angleyrad), 0, cos(angleyrad)); //y
		//rotationmatrix = (Mat_<float>(3,3) << 1, 0, 0, 0, cos(angleyrad), -sin(angleyrad), 0, sin(angleyrad), cos(angleyrad)); //x
		//rotationmatrix = (Mat_<float>(3,3) << cos(angleyrad), -sin(angleyrad), 0, sin(angleyrad), cos(angleyrad), 0, 0, 0, 1); //z
		
		for ( int i = 0; i < map_x.rows; i++ ) // here, i is for y and j is for x
			{
				for ( int j = 0; j < map_x.cols; j++ )
				{
					// normalizing to [-1, 1]
					xfish = (j - xcd) / halfcols;
					yfish = (i - ycd) / halfrows;
					rfish = sqrt(xfish*xfish + yfish*yfish);
					theta = atan2(yfish, xfish);
					phi = rfish*aperture/2;
					
					// Paul's co-ords - this is suitable when phi=0 is Pz=0
					
					//Px = cos(phi)*cos(theta);
					//Py = cos(phi)*sin(theta);
					//Pz = sin(phi);
					
					// standard co-ords - this is suitable when phi=pi/2 is Pz=0
					Px = sin(phi)*cos(theta);
					Py = sin(phi)*sin(theta);
					Pz = cos(phi);
					
					if(angley!=0 || anglex!=0)
					{
						// cos(angleyrad), 0, sin(angleyrad), 0, 1, 0, -sin(angleyrad), 0, cos(angleyrad));
						
						PxR = Px;
						PyR = cos(angleyrad) * Py - sin(angleyrad) * Pz;
						PzR = sin(angleyrad) * Py + cos(angleyrad) * Pz;
						
						Px = cos(anglexrad) * PxR - sin(anglexrad) * PyR;
						Py = sin(anglexrad) * PxR + cos(anglexrad) * PyR;
						Pz = PzR;
					}
					
					
					longi 	= atan2(Py, Px);
					lat	 	= atan2(Pz,sqrt(Px*Px + Py*Py));	
					// this gives south pole centred, ie yequi goes from [-1, 0]
					// Made into north pole centred by - (minus) in the final map_y assignment
					
					xequi = longi / CV_PI;
					// this maps to [-1, 1]
					yequi = 2*lat / CV_PI;
					// this maps to [-1, 0] for south pole
					
					if (rfish <= 1.1)		// outside that circle, let it be black
					// restored the black circle 
					// to help transformtype=5
					// avoid bottom pixels black, made it 1.1 instead of 1.0
					{
						map_x.at<float>(i, j) =  abs(xequi * map_x.cols / 2 + xcd);
						//map_y.at<float>(i, j) =  yequi * map_x.rows / 2 + ycd;
						// this gets south pole centred view
						
						// the abs is to correct for -0.5 xequi value at longi=0
						
						map_y.at<float>(i, j) =  yequi * map_x.rows / 2 + ycd;
						
					}
					
				 } // for j
				   
			} // for i
	 
     } // end of if transformtype == 0
	}	// end switch case
    
} // end function updatemap

std::string escaped(const std::string& input)
{
	// https://stackoverflow.com/questions/48260879/how-to-replace-with-in-c-string
    std::string output;
    output.reserve(input.size());
    for (const char c: input) {
        switch (c) {
            case '\a':  output += "\\a";        break;
            case '\b':  output += "\\b";        break;
            case '\f':  output += "\\f";        break;
            case '\n':  output += "\\n";        break;
            case '\r':  output += "\\r";        break;
            case '\t':  output += "\\t";        break;
            case '\v':  output += "\\v";        break;
            default:    output += c;            break;
        }
    }

    return output;
}

inline void writeIni(std::string iniwpath) 
{
	try {
		std::ofstream inifileout(iniwpath);
		inifileout << "#ini_file_for_OCVWarp--Comments_start_with_#" << std::endl;
		inifileout << "#Enter_each_parameter_in_the_line_below_the_comment. " << std::endl;
		inifileout << "#AngleXinDegrees_float" << std::endl;
		inifileout << anglexstr << std::endl;
		inifileout << "#AngleXIncrementperFrameinDegrees_float" << std::endl;
		inifileout << anglexincrstr << std::endl;
		inifileout << "#AngleYinDegrees_float" << std::endl;
		inifileout << angleystr << std::endl;
		inifileout << "#AngleYIncrementperFrameinDegrees_float" << std::endl;
		inifileout << angleyincrstr << std::endl;
		inifileout << "#Output_width_pixels" << std::endl;
		inifileout << outputw << std::endl;
		inifileout << "#Output_height_pixels" << std::endl;
		inifileout << outputh << std::endl;
		inifileout << "#0=360Fisheye_1=180Fisheye_etc--see_transformtype.txt" << std::endl;
		inifileout << transformtype << std::endl;
		inifileout << "#Output_video_codec_fourcc__use_NULL_for_same_as_input--see_fourcc.txt" << std::endl;
		inifileout << outputfourccstr << std::endl;		
		inifileout << "#Path_to_Map_file_used_for_transformtype_4_&_5" << std::endl;
		inifileout << "EP_xyuv_1920.map" << std::endl;
		inifileout << "#Output_fps_-1=same_as_input__0=image_sequence" << std::endl;
		inifileout << outputfpsstr << std::endl;
		
	} catch (int) {
		std::cerr << "An error occured writing to ini file."<< std::endl ; 		
	}
}


int main(int argc,char *argv[])
{
	////////////////////////////////////////////////////////////////////
	// Initializing variables
	////////////////////////////////////////////////////////////////////
	bool doneflag = 0, interactivemode = 0;
	bool showdisplay = 1;
    double anglex = 0;
    double angley = 0;
    double anglexincr = 0;
    double angleyincr = 0;
     
    int texturew = 2048;
    strpathtowarpfile = "EP_xyuv_1920.map";
    outputfourccstr[0] = 'N';
    outputfourccstr[1] = 'U';
    outputfourccstr[2] = 'L';
    outputfourccstr[3] = 'L';
    double outputfps = -1;  
    
    //const bool askOutputType = argv[3][0] =='Y';  // If false it will use the inputs codec type
    // this line above causes the windows build to not run! although it compiles ok.
    // askOutputType=1 works only on Windows (vfw?) currently
    const bool askOutputType = 0;

	// adding code to work with command-line arguments
	// argv[0] = name of exe, argv[1] = ini file path, argv[2] = input file path, argv[3] = output file path
	// but no error checking!
    char const * SaveFileName = "";
    char const * OpenFileNameini = "";
    char const * OpenFileName = "";
	bool argsSupplied = 0;
	if (argc == 4) {
		OpenFileNameini = argv[1];
		OpenFileName = argv[2];
		SaveFileName = argv[3];
		argsSupplied = 1;
		showdisplay = 0;
	}
		
    // adding code to open ini file instead of hardcoding
	// from https://github.com/hn-88/OCVvid2fulldome/
	std::string inistr;
	char const * FilterPatternsini[2] =  { "*.ini","*.*" };

		
	if(!argsSupplied) {	
	OpenFileNameini = tinyfd_openFileDialog(
				"Open an ini file if it exists",
				"",
				2,
				FilterPatternsini,
				NULL,
				0);
	}

	if (! OpenFileNameini) {
		// manual mode
		char const * lTmp;
		tinyfd_messageBox("Please Note", 
			"ini file not supplied or unreadable. So, manual inputs ...", 
			"ok", "info", 1);
		
		/* Data needed:
		 anglexstr;
		anglexincrstr;
		angleystr;
		angleyincrstr;
		outputw;
		outputh;
		transformtype;
		outputfourccstr;
		strpathtowarpfile;
                outputfpsstr
		*/
		lTmp = tinyfd_inputBox(
		"Please Input", "Output video width", "3840");
		if (!lTmp) return 1 ;	
		outputw = atoi(lTmp);
		
		lTmp = tinyfd_inputBox(
		"Please Input", "Output video height", "2160");
		if (!lTmp) return 1 ;	
		outputh = atoi(lTmp);

		// there are currently 3 input types and 4 output types
		// as seen in build/transformtype.txt

		// Input == Equirect --> Output can be 360fisheye=0, 180fisheye=1, warped=5
		// Input != Equirect --> Input can be 360fisheye or 180fisheye
		// Input == 360fisheye --> Output can be Equirect=2 only
		// Input == 180fisheye --> Output can be Equirect=3, warped=4
		
		int isInputFisheye180 = 0;
		int isOutputFisheye = 0;
		int isOutputFisheye180 = 0;
		int isOutputWarped = 0;
		
		int isInputEquirect = tinyfd_messageBox(
		"Transform type - Input" , /* NULL or "" */
		"Is the input Equirectangular VR360?"  , /* NULL or "" may contain \n \t */
		"yesno" , /* "ok" "okcancel" "yesno" "yesnocancel" */
		"question" , /* "info" "warning" "error" "question" */
		1 ) ;	/* 0 for cancel/no , 1 for ok/yes , 2 for no in yesnocancel */

		if (isInputEquirect != 1) {
			isInputFisheye180 = tinyfd_messageBox(
			"Transform type - Input" , 
			"Is the input 180 fisheye (fulldome)?"  , 
			"yesno" , 
			"question" , 
			1 ) ;	
		}

		isOutputFisheye = tinyfd_messageBox(
		"Transform type - Output" , 
		"Is the desired output fisheye?"  , 
		"yesno" , 
		"question" , 
		1 ) ;

		if (isOutputFisheye == 1) {
			isOutputFisheye180 = tinyfd_messageBox(
			"Transform type - Output" , 
			"Is the desired output 180 fisheye (fulldome)?"  , 
			"yesno" , 
			"question" , 
			1 ) ;	
		} else {
			isOutputWarped = tinyfd_messageBox(
			"Transform type - Output" , 
			"Is the desired output a warped file (mirrordome)?"  , 
			"yesno" , 
			"question" , 
			1 ) ;	
		}
		
		if (isInputEquirect == 1 ) {
			if (isOutputWarped == 1) transformtype = 5;
			else  {
			if ( isOutputFisheye180 == 1) transformtype = 1;
			else transformtype = 0;
			}
		} else if (isInputFisheye180 == 1) {
			if (isOutputWarped == 1) transformtype = 4;
			else  transformtype = 3;
		} else transformtype = 2;
		
		lTmp = tinyfd_inputBox(
		"Please Input", "Output FOURCC", "avc1");
		if (!lTmp) return 1 ;	
		std::strcpy(outputfourccstr,  lTmp);

		lTmp = tinyfd_inputBox(
		"Please Input", "Output fps (frames per second) - -1 to use input video fps, 0 for frame sequence", "-1");
		if (!lTmp) return 1 ;
		std::strcpy(outputfpsstr,  lTmp);
		
		// we need to ask user for anglex angley data only if isInputEquirect==1
		std::strcpy(anglexstr,"-90.0");
		std::strcpy(anglexincrstr,"0.0");
		std::strcpy(angleystr,"-160.0");
		std::strcpy(angleyincrstr,"0.0");
		std::strcpy(outputfpsstr,"-1.0");

		if (isInputEquirect==1) {
			lTmp = tinyfd_inputBox(
			"Please Input", "AngleX", "-90.0");
			if (!lTmp) return 1 ;	
			std::strcpy(anglexstr,  lTmp);
			
			lTmp = tinyfd_inputBox(
			"Please Input", "AngleX increment per frame", "0.0");
			if (!lTmp) return 1 ;	
			std::strcpy(anglexincrstr,  lTmp);
			
			lTmp = tinyfd_inputBox(
			"Please Input", "AngleY", "-160.0");
			if (!lTmp) return 1 ;	
			std::strcpy(angleystr,  lTmp);
			
			lTmp = tinyfd_inputBox(
			"Please Input", "AngleY increment per frame", "0.0");
			if (!lTmp) return 1 ;	
			std::strcpy(angleyincrstr,  lTmp);			
		}
		anglex = atof(anglexstr);
		angley = atof(angleystr);
		anglexincr = atof(anglexincrstr);
		angleyincr = atof(angleyincrstr);
		outputfps = atof(outputfpsstr);
	// here, we give an option for the user to save the ini file
        // if cancelled, the program just continues.
		char const * lIniFilterPatterns[1] = { "*.ini" };
    		char const * IniSaveFileName = tinyfd_saveFileDialog(
		"Choose the name and path of the ini file if you want to save the settings, like OCVWarp-4096ToW.ini",
		"OCVWarp-.ini",
		1,
		lIniFilterPatterns,
		NULL);
		if (IniSaveFileName)
		{
			writeIni(IniSaveFileName);
		}


	/* ***********
        adding a preview window would add significantly to complexity, so, skipping for now
		// Init cvui and tell it to create a OpenCV window, i.e. cv::namedWindow(WINDOW_NAME).
		cvui::init(WINDOW_NAME);
		cv::Mat frame = cv::Mat(cv::Size(400, 200), CV_8UC3);
		while (true) {
		// Fill the frame with a nice color
		frame = cv::Scalar(49, 52, 49);

		// Render UI components to the frame
		cvui::text(frame, 350, 10, "Preview");
	
		} // end while (true) loop
		************
         not adding a preview window for now
			*/
	} // if (! OpenFileNameini)
	
	else {
	     inistr = OpenFileNameini;
	     std::ifstream infile(inistr);
	    
	    
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
				infile >> anglexincrstr;
				infile >> tempstring;
				infile >> angleystr;
				infile >> tempstring;
				infile >> angleyincrstr;
				infile >> tempstring;
				infile >> outputw;
				infile >> tempstring;
				infile >> outputh;
				infile >> tempstring;
				infile >> transformtype;
				infile >> tempstring;
				infile >> outputfourccstr;
				infile >> tempstring;
				infile >> strpathtowarpfile;
				infile >> tempstring;
				infile >> outputfpsstr;
				infile.close();
				
				anglex = atof(anglexstr);
				angley = atof(angleystr);
				anglexincr = atof(anglexincrstr);
				angleyincr = atof(angleyincrstr);
				outputfps = atof(outputfpsstr);
			  }
	
		else std::cout << "Unable to open ini file, using defaults." << std::endl;
	
	} // end else block of if  (! OpenFileNameini) 
	
	std::cout << "Output codec type: " << outputfourccstr << std::endl;
	if(!argsSupplied) {
		namedWindow("Display", WINDOW_NORMAL | WINDOW_KEEPRATIO); // 0 = WINDOW_NORMAL
		resizeWindow("Display", round((float)(outputw)/(float)(outputh)*600), 600); // this doesn't work?
		moveWindow("Display", 0, 0);
	}
	
	char const * FilterPatterns[2] =  { "*.avi","*.*" };
	if(!argsSupplied) {
	 OpenFileName = tinyfd_openFileDialog(
		"Open a video file or image sequence with zero padded filenames",
		"",
		2,
		FilterPatterns,
		NULL,
		0);
	}

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
	// add image / image sequence support
	
	if (imread(OpenFileName).empty() ) {
		// https://docs.opencv.org/3.4/d4/da8/group__imgcodecs.html#ga288b8b3da0892bd651fce07b3bbd3a56
		// it is not an image file because imread returns empty, so must be a video file
	//std::cout  << "Input is not a single image. " << OpenFileName << std::endl;
	} // end if imread returns empty
	else {
		// check if it is an image sequence
		// https://stackoverflow.com/questions/538300/check-what-number-a-string-ends-with-in-c
		std::string test = OpenFileName;
		std::string::size_type pAt = test.find_last_of('.');                  // Find extension point
                std::string testwoext = test.substr(0, pAt);
		std::string::size_type last_char_pos = testwoext.find_last_not_of("0123456789");
		// std::cout  << "last not of is " << last_char_pos   << std::endl;
		std::string base = testwoext.substr(0, last_char_pos + 1);
		if (base == testwoext) {
			//std::cout  << "Input is not an image sequence. " << OpenFileName << std::endl;
			// Since input is a single image, 
			outputfps = 0;
		}else {
			// if it is an image sequence, change to the %0d format.
			// But don't offer to make any changes if commandline args are supplied.
			if(!argsSupplied) {
			        std::string::size_type last_num_pos = test.find_last_of("0123456789");
				// std::cout  << "last num pos is " << last_num_pos   << std::endl;
				// std::cout  << "Im seq " << base + "%0" + std::to_string(last_num_pos - last_char_pos) + "d" + test.substr(pAt)  << std::endl;
				// ask user before changing to %0d format
				bool isInputImSeq = tinyfd_messageBox(
				"Input type" , 
				"Is the input an image sequence with zero padded numbering starting from 0 or 1? If yes, we will change the path accordingly."  , 
				"yesno" , 
				"question" , 
				1 ) ;
				if(isInputImSeq) {
					std::string imgSeqPath = base + "%0" + std::to_string(last_num_pos - last_char_pos) + "d" + test.substr(pAt);
					OpenFileName = imgSeqPath.c_str();
				}
			}
		}
	}
	// reference:
	// https://docs.opencv.org/3.4/d7/d9e/tutorial_video_write.html
	
#ifdef __unix__
	
	VideoCapture inputVideo(OpenFileName);              // Open input
#else
// #ifdef _WIN64 assuming Windows
// otherwise, multiple definitions when building on Windows x64
	// Here, OpenCV on Windows needs escaped file paths. 
	// https://stackoverflow.com/questions/48260879/how-to-replace-with-in-c-string
	std::string escapedpath = escaped(std::string(OpenFileName));
	VideoCapture inputVideo(escapedpath.c_str());              // Open input
#endif
	
	if (!inputVideo.isOpened())
    {
        std::cout  << "Could not open the input video: " << OpenFileName << std::endl;
        return -1;
    }
     
#ifdef __unix__
    std::string OpenFileNamestr = OpenFileName; 
#else
//#ifdef _WIN64
	// Here, OpenCV on Windows needs escaped file paths. 
	std::string OpenFileNamestr = escapedpath; 
#endif 
	
  

    std::string::size_type pAt = OpenFileNamestr.find_last_of('.');                  // Find extension point
    std::string NAME = OpenFileNamestr.substr(0, pAt) + "F" + ".avi";   // Form the new name with container
    
    // here, we give an option for the user to choose the output file
    // path as well as type (container, like mp4, mov, avi).
	if(!argsSupplied) {
             SaveFileName = tinyfd_saveFileDialog(
		"Now enter the output video file name, like output.mp4",
		"",
		0,
		NULL,
		NULL);
	}

	if (! SaveFileName)
	{
		tinyfd_messageBox(
			"No output file chosen.",
			"Will be saved as inputfilename + F.avi",
			"ok",
			"info",
			1);
		 
	}
	else
	{
#ifdef __unix__
	NAME = std::string(SaveFileName);
#else
	// for Windows, escape the \ characters in the path
	escapedpath = escaped(std::string(SaveFileName));
	NAME = escapedpath;
#endif
	}
	
    
    int ex = static_cast<int>(inputVideo.get(CAP_PROP_FOURCC));     // Get Codec Type- Int form
    // Transform from int to char via Bitwise operators
    char EXT[] = {(char)(ex & 0XFF) , (char)((ex & 0XFF00) >> 8),(char)((ex & 0XFF0000) >> 16),(char)((ex & 0XFF000000) >> 24), 0};
    Size S = Size((int) inputVideo.get(CAP_PROP_FRAME_WIDTH),    // Acquire input size
                  (int) inputVideo.get(CAP_PROP_FRAME_HEIGHT));
    Size Sout = Size(outputw,outputh);            
    VideoWriter outputVideo;                                        // Open the output
//#ifdef _WIN64
	// OpenCV on Windows can ask for a suitable fourcc. 
	//outputVideo.open(NAME, -1, inputVideo.get(CAP_PROP_FPS), Sout, true);
	// this doesn't work well with the ffmpeg dll - don't use this.
//#endif 
		
    // if output fps is same as input fps, outputfps is set to -1,
	if (outputfps < 0) {
		outputfps = inputVideo.get(CAP_PROP_FPS);
	}
	
    if (!(outputfourccstr[0] == 'N' &&
    outputfourccstr[1] == 'U' &&
    outputfourccstr[2] == 'L' &&
    outputfourccstr[3] == 'L'))
        outputVideo.open(NAME, outputVideo.fourcc(outputfourccstr[0], outputfourccstr[1], outputfourccstr[2], outputfourccstr[3]), 
        outputfps, Sout, true);
    else
        outputVideo.open(NAME, ex, outputfps, Sout, true);


    if (!outputVideo.isOpened())
    {
        std::cout  << "Could not open the output video for write: " << NAME << std::endl;
        return -1;
    }
    std::cout << "Input frame resolution: Width=" << S.width << "  Height=" << S.height
         << " of nr#: " << inputVideo.get(CAP_PROP_FRAME_COUNT) << std::endl;
    std::cout << "Input codec type: " << EXT << std::endl;
     
    int  fps, key;
	int t_start, t_end;
    unsigned long long framenum = 0;
     
    Mat src, res, tmp;
    Mat dstfloat, dstmult, dstres, dstflip;
    
    std::vector<Mat> spl;
    Mat dst(Sout, CV_8UC3); // S = src.size, and src.type = CV_8UC3
    Mat dst2;	// temp dst, for double remap
    Mat map_x, map_y;
    if ((transformtype == 4) || (transformtype == 5) ) 
    {
		//~ map_x = Mat(Size(outputw*2,outputh*2), CV_32FC1);	// for 2x resampling
		//~ map_y = Mat(Size(outputw*2,outputh*2), CV_32FC1);
		// the above code causes gridlines to appear in output
		if (outputw<961)	//1K
			texturew = 1024;
			// debug
			//texturew = outputw;
		else if (outputw<1921)	//2K
			texturew = 2048;
		else if (outputw<3841)	//4K
			texturew = 4096;
		else // (outputw<7681)	//8K
			texturew = 8192;
			// debug - had set Size to outputw,h
		map_x = Mat(Size(texturew,texturew), CV_32FC1);	// for upsampling
		map_y = Mat(Size(texturew,texturew), CV_32FC1);
		if (transformtype == 5)
		{
			dst2 = Mat(Size(texturew,texturew), CV_32FC1);
		}
	}
	else
	{
		map_x = Mat(Sout, CV_32FC1);
		map_y = Mat(Sout, CV_32FC1);
	}
    Mat dst_x, dst_y;
    // Mat map2x, map2y // these are made global vars 
    Mat dst2x, dst2y;	// for transformtype=5, double remap
    if (transformtype == 5) 
    {
		map2x = Mat(Size(texturew,texturew), CV_32FC1);	
		map2y = Mat(Size(texturew,texturew), CV_32FC1);
		map2x = Scalar((texturew+texturew)*10);
		map2y = Scalar((texturew+texturew)*10);
		// initializing so that it points outside the image
		// so that unavailable pixels will be black
	}
	
    map_x = Scalar((outputw+outputh)*10);
    map_y = Scalar((outputw+outputh)*10);
    // initializing so that it points outside the image
    // so that unavailable pixels will be black
    
    if (transformtype == 5) 
    {
		update_map(anglex, angley, map_x, map_y, 4);
		// first initialize map_x and map_y for final warp,
		// which is exactly like transformtype=4
		update_map(anglex, angley, map2x, map2y, 5);
		convertMaps(map_x, map_y, dst_x, dst_y, CV_16SC2);	// supposed to make it faster to remap
		convertMaps(map2x, map2y, dst2x, dst2y, CV_16SC2);
		
	}
	else
	{
		update_map(anglex, angley, map_x, map_y, transformtype);
		// debug
		//~ dst_x = map_x;
		//~ dst_y = map_y;
		convertMaps(map_x, map_y, dst_x, dst_y, CV_16SC2);	// supposed to make it faster to remap
	}
    
    t_start = time(NULL);
	fps = 0;
	
    for(;;)
    {
        inputVideo >> src;              // read
        if (src.empty()) break;         // check if at end
        //imshow("Display",src);
        key = waitKey(10);
        
        if(interactivemode)
        {
			if (transformtype==5)
			{
				// only map2 needs to be updated
			update_map(anglex, angley, map2x, map2y, 5);
			convertMaps(map2x, map2y, dst2x, dst2y, CV_16SC2);	// supposed to make it faster to remap
			}
			else
			{
			update_map(anglex, angley, map_x, map_y, transformtype);
			convertMaps(map_x, map_y, dst_x, dst_y, CV_16SC2);	// supposed to make it faster to remap
			}
    
			interactivemode = 0;
		}
		
		switch (transformtype)
				{
				case 5: // 360 to 180 fisheye and then to warped
					resize( src, res, Size(texturew, texturew), 0, 0, INTER_CUBIC);
					break;

				case 4: // 180 fisheye to warped
					// the transform needs a flipped source image, flipud
					flip(src, src, 0);	// because the mesh assumes 0,0 is bottom left
					//debug - had changed to outputw, h
					resize( src, res, Size(texturew, texturew), 0, 0, INTER_CUBIC);
					break;

				case 3: // 360 fisheye to Equirect
					resize( src, res, Size(outputw, outputh), 0, 0, INTER_CUBIC);
					break;
					
				case 2: // 360 fisheye to Equirect
					resize( src, res, Size(outputw, outputh), 0, 0, INTER_CUBIC);
					break;
					
				case 1: // Equirect to 180 fisheye
					resize( src, res, Size(outputw, outputh), 0, 0, INTER_CUBIC);
					break;
				
				default:	
				case 0: // Equirect to 360 fisheye
					resize( src, res, Size(outputw, outputh), 0, 0, INTER_CUBIC);
					break;
					
				}
		if (transformtype == 5)
		{
			// here we have two remaps
			remap( res, dst2, dst2x, dst2y, INTER_LINEAR, BORDER_CONSTANT, Scalar(0, 0, 0) );
			// the second remap needs flipping, like transformtype=4
			flip(dst2, dst2, 0);	// because the mesh assumes 0,0 is bottom left
			remap( dst2, dst, dst_x, dst_y, INTER_LINEAR, BORDER_CONSTANT, Scalar(0, 0, 0) );
		}
		else
		{
        remap( res, dst, dst_x, dst_y, INTER_LINEAR, BORDER_CONSTANT, Scalar(0, 0, 0) );
		}
	
        if ((transformtype == 4) || (transformtype == 5) )
        {
			// multiply by the intensity Mat
			dst.convertTo(dstfloat, CV_32FC3);
			multiply(dstfloat, I, dstmult);
			//debug
			//dstmult = dstfloat;
			dstmult.convertTo(dstres, CV_8UC3);
			// this transform is 2x2 oversampled
			resize(dstres, dstflip, Size(outputw,outputh), 0, 0, INTER_AREA);
			flip(dstflip, dst, 0); 	// flip up down again
		}
			
        if(showdisplay)
			imshow("Display", dst);
			
        printf("\r");
        
        fps++;
        t_end = time(NULL);
		if (t_end - t_start >= 5)
		{
			printf("Frame: %llu x: %.0f y: %.0f fps: %.1f           \r", framenum++, anglex, angley, float(fps)/5 );
			// extra spaces to delete previous line's characters if any
			fflush(stdout);
			t_start = time(NULL);
			fps = 0;
		}
		else
		{
			printf("Frame: %llu x: %.0f y: %.0f \r", framenum++, anglex, angley );
			fflush(stdout);
		}
			
        
       //outputVideo.write(res); //save or
       outputVideo << dst;
       if (anglexincr!=0 || angleyincr!=0) {
	   anglex = anglex + anglexincr;
	   angley = angley + angleyincr;
	   interactivemode = 1;
       }
	
       
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
					interactivemode = 1;
					break;
					
				case 'm':
				case '-':
				case '_':	// decrease angley
					angley = angley - 1.0;
					interactivemode = 1;
					break;
					
				case 'k':
				case '}':
				case ']':	// increase anglex
					anglex = anglex + 1.0;
					interactivemode = 1;
					break;
					
				case 'h':
				case '{':
				case '[':	// decrease anglex
					anglex = anglex - 1.0;
					interactivemode = 1;
					break;
				
				case 'U':
					// increase angley
					angley = angley + 10.0;
					interactivemode = 1;
					break;
					
				case 'M':
					// decrease angley
					angley = angley - 10.0;
					interactivemode = 1;
					break;
					
				case 'K':
					// increase anglex
					anglex = anglex + 10.0;
					interactivemode = 1;
					break;
					
				case 'H':
					// decrease anglex
					anglex = anglex - 10.0;
					interactivemode = 1;
					break;	
					
				case 'D':
				case 'd':
					// toggle display
					if(showdisplay)
						showdisplay=0;
					else
						showdisplay=1;
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
