#include "header.h"

using namespace std;


extern long seed;
extern long int tim;

extern struct fibers parameters[numberOfDetectors+1];

void get_media_morphology();

double random(long *idum);

//Initialize position and direction of light with consideration of NA and angle of fiber
//=========================================================================================================================================================================================================================================
void Launch(double &x,double &y,double &z,double &ux,double &uy,double &uz) //初始化光子一開始的位置、方向、權重 
{
	 double rnd_position;

	 double rnd_Azimuth;

	 double rnd_direction;

	 double rnd_rotated;

	 double AzimuthAngle;

	 double launchPosition;

	 double theta_direction;

	 double rotated_angle;

	 get_media_morphology();

	 rnd_position=rand()/(double)RAND_MAX;

     rnd_Azimuth=rand()/(double)RAND_MAX;
 
	 rnd_direction=rand()/(double)RAND_MAX;

	 rnd_rotated=rand()/(double)RAND_MAX;

	 AzimuthAngle=2*PI*rnd_Azimuth;

	 rotated_angle=2*PI*rnd_rotated;


     launchPosition=parameters[0].radius*sqrt(rnd_position);

	 x=launchPosition*cos(AzimuthAngle)/cos(parameters[0].angle);

     y=launchPosition*sin(AzimuthAngle);
     
     z=0.0;

	 theta_direction=asin(parameters[0].NA/n_source)*rnd_direction;

	 
     ux=(sin(theta_direction)*cos(rotated_angle)*cos(parameters[0].angle))
		-(cos(theta_direction)*sin(parameters[0].angle));

     uy=sin(theta_direction)*sin(rotated_angle);

     uz=(sin(theta_direction)*cos(rotated_angle)*sin(parameters[0].angle))
		+(cos(theta_direction)*cos(parameters[0].angle));

	 if(detect==0)
	 {
		 x=0,y=0,z=0;
	 }


}
//=========================================================================================================================================================================================================================================
