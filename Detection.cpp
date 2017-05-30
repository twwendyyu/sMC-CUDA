#include "header.h"
using namespace std;


extern struct fibers parameters[numberOfDetectors+1];
extern struct ArrayOfPhotonPath photon_record;


//check for direction of photon whether smaller than NA of Detector or not
bool DetectionOfPhoton(double &ux, double &uz)
{
	int detected=numberOfDetectors+1;

	bool detector_identified=false;	

	double uz_rotated;

	for (int i=1; i<=numberOfDetectors && !detector_identified; i++)
	{
		uz_rotated=(ux*sin(parameters[i].angle))+(uz*cos(parameters[i].angle)); 
		
		if (acos(fabs(uz_rotated))<=asin(parameters[i].NA/ n_detector))
		{
			
			detector_identified=true;
		}	
	}

	return detector_identified;
}
