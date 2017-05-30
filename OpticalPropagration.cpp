#include <iostream>
#include <stdio.h>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>
#include "header.h"

using namespace std;

struct ArrayOfPhotonPath photon_record;

extern struct fibers parameters[numberOfDetectors+1];

extern struct media morphology[numberOfMedium+2];

extern long seed;
extern long int tim;
double random(long *idum);

void RecordPath(double &,double & ,double &);

bool  Specular(double &uz);
void ModifyPositionAndDirection(double &x,double &y,double &z,double &ux,double &uy,double &uz,int &layer,bool TIF);
bool checkForTIF(double &x,double &y,double &z,double &ux,double &uy,double &uz,int &layer,double &weight,double &WR, double &WA,double &WT);
void changeDirectionDueToRefraction(double &ux,double &uy,double &uz,int &layer);
double fresnel(double &n_in,double &n_out,double &uz);
void InitializationOfWeight(double &weight);
void getUnitStep(double & step_length,int &layer) ;
void move(double & x,double & y,double & z , double & ux,double & uy,double & uz,double & stepLength);
void reflect(double &WR,double &weight)   ;
void transmit(double &WT,double &weight) ;
int sign(double input) ;
void direction(double &ux,double &uy,double &uz,int &layer);
bool roulette(double &z,double &ux,double &uy,double &uz,double &weight,int &layer);
void absorption(double &WA,int &layer,double &weight);
void newPosition(double &z,double &uz,int &layer);
void ClearRecordPath();
//============================================================================================================
void  get_fiber_parameters();
void  get_media_morphology();
//============================================================================================================
bool  Specular(double &uz)
{

	bool specular_reflection=false;	/* Default value is 0; If internal reflection will	occur, becomes 1 */

	double theta_c;		/* Critical angle for incidence */
	
	double rnd;		/* Fresnel reflectance */

	double r;

	double n_in=0,n_out=0;
	
	n_in=morphology[0].n;
	n_out=morphology[1].n;
	
	if (n_out<n_in)	// Total internal reflection possible
	{
		theta_c=asin(n_out/n_in);
	}
	else
	{
		theta_c=PI/2.0;

	}

	specular_reflection=false;
        
	if (acos(fabs(uz))>=theta_c)	// Total internal reflection       
	{
		specular_reflection=true; 
	}   
	else if (acos(fabs(uz))<theta_c)
    {
		r=fresnel(n_in,n_out,uz );
		rnd=rand()/(double)RAND_MAX;	/* Generates a random number in (0, 1) */
        if (rnd<=r)			// Internal reflection
        {
			specular_reflection=true;
		}
	}
	return specular_reflection;

}
void ModifyPositionAndDirection(double &x,double &y,double &z,double &ux,double &uy,double &uz,int &layer,bool TIF)
{
	double DIB=0;

	if(uz<0)
	{
		DIB=(0-z)/uz;
	
	}
	else if(uz>0)
	{
		DIB=(morphology[1].Thickness-z)/uz;
		
	}
	if(TIF)
	{
		photon_record.Array_site-=1;
		
		move(x,y,z,ux,uy,uz,DIB);
		
		DIB=-DIB;

		uz=-uz;
	
		move(x,y,z,ux,uy,uz,DIB);
	
	}
	else
	{

		photon_record.Array_site-=1;

		move(x,y,z,ux,uy,uz,DIB);

		changeDirectionDueToRefraction(ux,uy,uz,layer);
	}

}
void RecordPath(double &x,double &y ,double &z)
{
	  photon_record.Array[0][photon_record.Array_site]=x;
	  photon_record.Array[1][photon_record.Array_site]=y;
	  photon_record.Array[2][photon_record.Array_site]=z;

      photon_record.Array_site++;


}
void ClearRecordPath() 
{
	 photon_record.Array_site-=photon_record.Array_site;
}

bool checkForTIF(double &x,double &y,double &z,double &ux,double &uy,double &uz,int &layer,double &weight,double &WR, double &WA,double &WT)
{
	int photonOut=0;	/* Default value is 0; If internal reflection will	occur, becomes 1 */

	double theta_c;		/* Critical angle for incidence */
	
	double rnd;		/* Fresnel reflectance */

	double r;
	double n_in=0,n_out=0;
	bool TIF=false;
	if(uz>0)
	{
		n_in=morphology[1].n;
		n_out=morphology[2].n;
	}
	else if(uz<0)
	{
		n_in=morphology[1].n;
		n_out=morphology[0].n;
	}
	/* Calculates the critical angle */
	if (n_out<n_in)	// Total internal reflection possible
	{
        
		theta_c=asin(n_out/n_in);
        if (acos(fabs(uz))>=theta_c)	// Total internal reflection
        {
			TIF=true;
        }
		else
		{
			 r=fresnel(n_in,n_out,uz );
             rnd=rand()/(double)RAND_MAX;	/* Generates a random number in (0, 1) */
             if (rnd<=r)			// Internal reflection
             {
				 TIF=true;
             }
		}
    }  
    else 
    {            
             r=fresnel(n_in,n_out,uz );
             rnd=rand()/(double)RAND_MAX;	/* Generates a random number in (0, 1) */
             if (rnd<=r)			// Internal reflection
             {
				 TIF=true;
             }
    }
	return TIF;

}


void changeDirectionDueToRefraction(double &ux,double &uy,double &uz,int &layer)
{
	
	double theta_i, theta_t;	/* Incident and transmitted angles, respectively */
	
	double n_in=0,n_out=0;

	if(uz>0)
	{
		n_in=morphology[0].n;
		n_out=morphology[1].n;
	}
	else if(uz<=0)
	{
		n_in=morphology[1].n;
		n_out=morphology[0].n;
	}

	theta_i=acos( fabs(uz) );
	theta_t=asin( (n_in*sin(theta_i))/n_out);
	
	ux=(ux)*(n_in/n_out);
	uy=(uy)*(n_in/n_out);
	uz=(sign(uz))*cos(theta_t);

}
void InitializationOfWeight(double &weight)
{
	 double R_sp; 

     R_sp=((n_source-morphology[1].n)*(n_source-morphology[1].n))/((n_source+morphology[1].n)*(n_source+morphology[1].n));


     weight=1-R_sp;
     
}
void getUnitStep(double & step_length,int &layer)     
{

     double rn=0.0;

	 do
	 { 
		rn=rand()/(double)RAND_MAX;

	 }while(rn==0);

	 step_length=(-1*log(rn)/(double)(morphology[1].Absorb_of_Baseline+morphology[1].Scattering_of_Baseline));
		
}

void move(double & x,double & y,double & z , double & ux,double & uy,double & uz,double & stepLength)
{
	x=x+stepLength*ux;
    y=y+stepLength*uy;
	z=z+stepLength*uz;

	RecordPath(x,y,z);

	return;
}


void reflect(double &WR,double &weight)   
{
	
     WR=WR+weight;
     return;
}

void transmit(double &WT,double &weight) 
{ 
     WT=WT+weight;
     return;
}


int sign(double input) 
{

	if (input>0)
	{
		return 1;
	}
	
	else if (input<0)
	{
		return -1;
	}

	else
	{
		return 0;
	}

} 
void direction(double &ux,double &uy,double &uz,int &layer)  
{     
      double rn,temp,theta,phi;

	  double uxprime,uyprime,uzprime;
      
      double g=0;

	   g=morphology[1].Anisotropy;

       if (g==0)
       {
		   rn=rand()/(double)RAND_MAX;
		   theta = acos(2*rn-1);
	   }
	   else
       { 
		   rn=rand()/(double)RAND_MAX;
		   temp = (1-g*g)/(1-g +2*g*rn);
		   theta=acos((1/(2*g))*(1 + g*g - (temp*temp)));
	   }
        
       rn=rand()/(double)RAND_MAX;
       
       phi = 2*PI*rn;
       
       if(abs(uz)>0.9999)
	   {
		   uxprime=sin(theta)*cos(phi);
		   uyprime=sin(theta)*sin(phi);
		   uzprime=uz*cos(theta)/abs(uz);
	   }
       else
	   {
		   uxprime= sin(theta)/sqrt(1-uz*uz)*(ux*uz*cos(phi)-uy*sin(phi)) + ux*cos(theta); 
		   uyprime= sin(theta)/sqrt(1-uz*uz)*(uy*uz*cos(phi)+ux*sin(phi)) + uy*cos(theta);
		   uzprime= -sin(theta)*cos(phi)*sqrt(1-((uz)*(uz))) + uz*cos(theta);
	   } 

	   ux=uxprime;
       uy=uyprime;
       uz=uzprime;

	   return ;

}
bool roulette(double &z,double &ux,double &uy,double &uz,double &weight,int &layer)       //判斷 光子是否存活  (當光子被吸收到很小時才會需要做此判斷 (在test_in_tissue會用到)) 
{
      double rn=rand()/(double)RAND_MAX;   
	  bool absorb=false;
        if (rn<(1/m))
        { weight=m*weight;
          direction(ux,uy,uz,layer);
        }
        else
        {weight=0;absorb=true;}

		return absorb;
}

void absorption(double &WA,int &layer,double &weight)
{
	WA=WA+weight*(morphology[layer].Absorb_of_Baseline/(morphology[layer].Absorb_of_Baseline+morphology[layer].Scattering_of_Baseline));

	weight=weight*(morphology[layer].Scattering_of_Baseline/(morphology[layer].Absorb_of_Baseline+morphology[layer].Scattering_of_Baseline));
}

void newPosition(double &z,double &uz,int &layer)
{
     if (z<0)
     {z=-z;uz=-uz;}
     else if(z >morphology[1].Thickness)
	 {z=2*morphology[1].Thickness-z;uz=-uz;}
}
double fresnel(double &n_in,double &n_out,double &uz)                             
{    
     double theta_t,R,theta_in;

     theta_in=acos(fabs(uz));

     theta_t = asin((n_in/n_out)*sin(theta_in));

     R=0.5 * ( (sin(theta_in-theta_t)*sin(theta_in-theta_t)) / (sin(theta_in+theta_t)*sin(theta_in+theta_t))+(tan(theta_in-theta_t)*tan(theta_in-theta_t)) / (tan(theta_in+theta_t)*tan(theta_in+theta_t)) );

     return R;
}




     
