#include "header.h"
using namespace std;

struct fibers parameters[numberOfDetectors+1];

//Set the parameters of source fiber and detector

void get_fiber_parameters()
{
		
	parameters[0].radius=illumination_r	;			
	parameters[0].NA=NAOfSource;				
	parameters[0].position=0.0;			
	parameters[0].angle=AngleOfSource*PI/180;
	
	for(int Number=1;Number<=numberOfDetectors;Number++)
	{
            parameters[Number].radius=collect_r;	

            parameters[Number].NA=NAOfDetector;		
		
			parameters[Number].angle=AngleOfDetector*PI/180;

			parameters[Number].position=((illumination_r)+(collect_r)*(2*Number-1))/cos(AngleOfDetector*PI/180);	
    }
	
	
	return ;

}      

//illumination_r+collect_r*(2*Number-1)
 
