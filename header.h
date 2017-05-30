#include <cmath>
#include <vector>
#include <string>
#include <sstream> 
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

#define NumberOfPhotons 10000000

#define TotalThickness 1.0

#define Threshold 2000


//在此區調整光纖參數

#define n_source 1.61

#define n_detector 1.61

#define n_Medium 1.4

#define NAOfSource 0.4

#define NAOfDetector 0.4

#define AngleOfSource -46.0 //toward x axis :positive

#define AngleOfDetector 46.0 //toward x axis :negative

#define detect 1// convolution:0 another:1





#define PI 3.14159265

#define m 20

#define w_th 0.005

#define numberOfDetectors 10

#define numberOfMedium 2

#define LamdaNumber  376

#define illumination_r 0.01

#define collect_r 0.01

//程式應用結構
struct ArrayOfPhotonPath 
{
       int Array_site;

       double Array[3][Threshold];

	   double vector[2];
};

struct GridSystemForReflectance
{
	double Array[100][100];  
};

struct GridSystemForAbsorbance
{
	double Array[100][100];

};

struct ScalingResult
{        
	double weight[LamdaNumber];
	double r[LamdaNumber];
	double r_vs_weight[LamdaNumber][10];
};

struct mua_array
{double all[LamdaNumber];};

struct mus_array
{double all[LamdaNumber];};

struct crash
{double crash1,crash2;};

struct photon
{double x,y,z,weight;};  

struct media 
{
	double Thickness ;
	double n;
    double Anisotropy;
    double Absorb_of_Baseline;
    double Scattering_of_Baseline;        
};

struct fibers 
{
	double radius;			/* Radius of the fiber in cm */
	double NA;				/* NA of the fiber (in air) */
	double position;		/* Distance from the origin to the center position of the 
							fiber in cm (defined to be 0.0 for the source fiber) */			
	double angle;			/* Tilt angle of the fiber in radians */
};
