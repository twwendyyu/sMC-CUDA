#include "header.h"

struct media morphology[numberOfMedium+2];


//Set up morphology of each layer
void get_media_morphology()
{
	morphology[0].Thickness=0.0;				
	morphology[0].n=n_source;
	morphology[0].Absorb_of_Baseline=0.0;
	morphology[0].Scattering_of_Baseline=0.0;
	morphology[0].Anisotropy=0.9;
   
	
	morphology[1].Thickness=TotalThickness;				
	morphology[1].n=n_Medium;
	morphology[1].Absorb_of_Baseline=0.0;
	morphology[1].Scattering_of_Baseline=200.0;
	morphology[1].Anisotropy=0.9;   

	morphology[2].Thickness=0.0;				
	morphology[2].n=n_Medium;
	morphology[2].Absorb_of_Baseline=0.0;
	morphology[2].Scattering_of_Baseline=0.0;
	morphology[2].Anisotropy=0.9;
}





