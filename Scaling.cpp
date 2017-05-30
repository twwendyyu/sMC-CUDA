#include <iostream>
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include <math.h>
#include <fstream>
#include <string>
#include <float.h>
#include "header.h"

#define mut_0 (morphology[1].Absorb_of_Baseline+morphology[1].Scattering_of_Baseline)
#define afa_0 (morphology[1].Scattering_of_Baseline/(morphology[1].Absorb_of_Baseline+morphology[1].Scattering_of_Baseline))
#define mut_1 (ms_1+ma_1)
#define afa_1 (ms_1/mut_1)
#define mut_2 (ms_2+ma_2)
#define afa_2 (ms_2/mut_2)
using namespace std;

extern struct fibers parameters[numberOfDetectors+1];
extern struct media morphology[numberOfMedium+2];
extern struct ArrayOfPhotonPath photon_record;
extern struct ScalingResult data;
extern struct ScalingResult data2;
extern struct ScalingResult data3;

void scale(double ma_1,double ms_1,double ma_2,double ms_2,double thickness1,struct photon &_new,int &line_number,struct photon &_old,struct photon &_out,struct crash &_count);
void shift(double layer,struct photon &old,struct photon &_new,struct photon &_out,struct crash &c,double ma_1,double ms_1,double ma_2,double ms_2);
double p(double x,double Collect2Source,double r);
double simpson(double a, double b, int n,double Collect2Source,double r);

void  get_fiber_parameters();
void  get_media_morphology();

void Scaling(struct  mua_array &ma1_out,struct  mus_array &ms1_out,struct  mua_array &ma2_out,struct  mus_array &ms2_out,float &Thickness_UpLayer)
{
	
    int as=0,as1=0,as2=0,as3=0,as4=0,as5=0,as6=0,as7=0,as8=0,as9=0;
    int bs=0,bs1=0,bs2=0,bs3=0,bs4=0,bs5=0,bs6=0,bs7=0,bs8=0,bs9=0;

    double weight,ma_1=0.0,ma_2=0.0,ms_1=0.0,ms_2=0.0;
	int detector_identified=0;
	int NeedRecord=0;
    int line_number=0;
    double Collect2Source=0.0;
    double lo=0.0,up=0.0,r=0.0;  
    double p[10];
	for(int pI=0;pI<10;pI++)
	{
		p[pI]=0;
	}
	int CounterOfPhotons=0;
	int run_number=0;
	bool Flag=false;
	get_fiber_parameters();
	get_media_morphology();

	photon _new;
    photon old;
    photon out;
	crash _count;
	ma_1=0;
	ms_1=0;
	ma_2=0;
	ms_2=0;
	_new.x=0, _new.y=0, _new.z=0;	
	old.x=0, old.y=0, old.z=0;  
	out.x=0, out.y=0, out.z=0;
	/*ofstream cc;
	cc.open("cc.txt");*/

	for(int i=75;i<276;i++) 
	{
		     detector_identified=0;
			 _count.crash1=0;
			 _count.crash2=0;
			 line_number=0;			 
			 CounterOfPhotons=0;
			 run_number= photon_record.Array_site;
			 
			 ma_1 = ma1_out.all[i];
			 ms_1 = ms1_out.all[i];
			 ma_2 = ma2_out.all[i];
			 ms_2 = ms2_out.all[i];


			 Flag=true;
			 while(Flag)
			 {
				 if(run_number==CounterOfPhotons)          
				 {Flag=false;}
				 else
				 { 
					 _new.x=photon_record.Array[0][CounterOfPhotons];
					 _new.y=photon_record.Array[1][CounterOfPhotons];
					 _new.z=photon_record.Array[2][CounterOfPhotons];

					 // cc<<setw(14)<<_new.x<<setw(14)<<_new.y<<setw(14)<<_new.z<<"   "<<run_number<<endl;

					 scale(ma_1,ms_1,ma_2,ms_2,Thickness_UpLayer,_new,line_number,old,out,_count);

                     CounterOfPhotons++;  

					 Flag=true;    
				 } 
              } 
			
			 weight=1*(pow(afa_1/afa_0,_count.crash1-1))*(pow(afa_2/afa_0,_count.crash2));

			 //cout<<out.z<<endl;

			 r=sqrt(pow(out.x,2)+pow(out.y,2)); 

			 data.r[i]=r;    

			 data.weight[i]+=weight; 

			 double FirstItem=0.0,SecondItem=0.0,ThirdItem=0.0;

			 if(detect==0)
			 {
				 for(int j=0;j<numberOfDetectors;j++)        
				 {
					 /*if(r>(j+1)*0.01 && r<(j+2)*0.01)
					 data.r_vs_weight[i][j]+=weight;
				 }*/

                   /* Collect2Source=parameters[j+1].position;
                    lo=max(-illumination_r,Collect2Source-r-collect_r);
                    up=min(illumination_r,Collect2Source-r+collect_r);                  
                    if (Collect2Source>=(illumination_r+collect_r)&& (r>=(Collect2Source-illumination_r-collect_r))&&(r<(Collect2Source+illumination_r+collect_r)))
                    {p[j]=(2/(PI*PI*illumination_r*illumination_r))*simpson(lo,up,10,Collect2Source,r);}
                    else 
                    {p[j]=0;}
                    if(p[j]<= DBL_MAX && p[j] >= -DBL_MAX)
                    {
						data.r_vs_weight[i][j]+=(weight*p[j]);
					}
                    else
                    {data.r_vs_weight[i][j]+=0;}  
				 }*/
					p[j]=0.0;
                    Collect2Source=parameters[j+1].position;
                    lo=max(-illumination_r,Collect2Source-r-collect_r);
                    up=min(illumination_r,Collect2Source-r+collect_r);                  
                    if (Collect2Source>=(illumination_r+collect_r)&& (r>=(Collect2Source-illumination_r-collect_r))&&(r<(Collect2Source+illumination_r+collect_r)))
                    {
						for(double x=lo;x<up;x+=0.00005)
						{
							double FirstItem=0.0,SecondItem=0.0,ThirdItem=0.0;
							FirstItem=(Collect2Source-x);

							SecondItem=acos((Collect2Source*Collect2Source+(Collect2Source-x)*(Collect2Source-x)-illumination_r*illumination_r)/(2*(Collect2Source-x)*Collect2Source)) ;
							if((Collect2Source*Collect2Source+(Collect2Source-x)*(Collect2Source-x)-illumination_r*illumination_r)/(2*(Collect2Source-x)*Collect2Source)>=1)
							{SecondItem=0.0;}

							ThirdItem=acos((r*r+(Collect2Source-x)*(Collect2Source-x)-collect_r*collect_r)/(2*(Collect2Source-x)*r));

							if((r*r+(Collect2Source-x)*(Collect2Source-x)-collect_r*collect_r)/(2*(Collect2Source-x)*r)>=1)
							{ThirdItem=0.0;}
					        
							p[j]+=FirstItem*SecondItem*ThirdItem;
						}
						p[j]*=(2/(PI*PI*illumination_r*illumination_r))*0.00005;
					}
                    else 
                    {p[j]=0;}
                    if(p[j]<= DBL_MAX && p[j] >= -DBL_MAX)
                    {
						data.r_vs_weight[i][j]+=(weight*p[j]);
					}
                    else
                    {
						data.r_vs_weight[i][j]+=0;
					} 
				 }
			 }

			 else if(detect==1)
			 {
				for (int f=1; f<=numberOfDetectors ; f++)
				{
					if ( ( ((out.x-parameters[f].position)/(parameters[f].radius/cos(parameters[f].angle)))*
						((out.x-parameters[f].position)/(parameters[f].radius/cos(parameters[f].angle))) ) + 
						( ((out.y)/parameters[f].radius)*((out.y)/parameters[f].radius) ) <=1.0)
					{
						data.r_vs_weight[i][f-1]+=weight;
						detector_identified=1;

						/*for(int q=0;q<photon_record.Array_site;q++)
						{
							for(int qq=0;qq<3;qq++)
							{
								aa<<setw(14)<<photon_record.Array[qq][q];
							}
							aa<<endl;
						}
						aa<<endl;*/
					}
				}

			 }
	}
	   
	   return ;

}
void shift(double layer,struct photon &old,struct photon &_new,struct photon &_out,struct crash &c,double ma_1,double ms_1,double ma_2,double ms_2)
{   
     if (old.z<=layer)
     {
		 c.crash1+=1;          
                    if (_new.z<=layer)
                    {
						            
                                     _out.x=_out.x+(mut_0/mut_1)*(_new.x-old.x);
                                     _out.y=_out.y+(mut_0/mut_1)*(_new.y-old.y);
                                     _out.z=_out.z+(mut_0/mut_1)*(_new.z-old.z);
                    }
                    else 
                    {                //定義在peusudo layer的位置 
                                     photon p;
                                     p.x=old.x+(layer-old.z)/(_new.z-old.z)*(_new.x-old.x);
                                     p.y=old.y+(layer-old.z)/(_new.z-old.z)*(_new.y-old.y);
                                     p.z=layer;
                                     //===================================================================
                                     //step3 
                                     //縮放在每一層pseudo layer內的水平位移，加總起來得到總共的水平位移
                                     _out.x=_out.x+(mut_0/mut_1)*(p.x-old.x);
                                     _out.y=_out.y+(mut_0/mut_1)*(p.y-old.y);
                                     _out.z=_out.z+(mut_0/mut_1)*(p.z-old.z);
                                     
                                     _out.x=_out.x+(mut_0/mut_2)*(_new.x-p.x);
                                     _out.y=_out.y+(mut_0/mut_2)*(_new.y-p.y);
                                     _out.z=_out.z+(mut_0/mut_2)*(_new.z-p.z);
							
                                     //===================================================================
					}
     }
     else 
     {
		 c.crash2+=1;   
                    if (_new.z>=layer)
                    {
                                     _out.x=_out.x+(mut_0/mut_2)*(_new.x-old.x);
                                     _out.y=_out.y+(mut_0/mut_2)*(_new.y-old.y);
                                     _out.z=_out.z+(mut_0/mut_2)*(_new.z-old.z);
                    }
                    
                    else 
                    {                //定義在peusudo layer的位置 
                                     photon p;
                                     p.x=old.x+(layer-old.z)/(_new.z-old.z)*(_new.x-old.x);
                                     p.y=old.y+(layer-old.z)/(_new.z-old.z)*(_new.y-old.y);
                                     p.z=layer;
                                     //=================================================================
                                     //step3 
                                     //縮放在每一層pseudo layer內的水平位移，加總起來得到總共的水平位移
                                     _out.x=_out.x+(mut_0/mut_2)*(p.x-old.x);
                                     _out.y=_out.y+(mut_0/mut_2)*(p.y-old.y);
                                     _out.z=_out.z+(mut_0/mut_2)*(p.z-old.z);
                                     
                                     _out.x=_out.x+(mut_0/mut_1)*(_new.x-p.x);
                                     _out.y=_out.y+(mut_0/mut_1)*(_new.y-p.y);
                                     _out.z=_out.z+(mut_0/mut_1)*(_new.z-p.z);
                                     //==================================================================
									
                                                   
                                     }
            }
}

void scale(double ma_1,double ms_1,double ma_2,double ms_2,double thickness1,struct photon &_new,int &line_number,struct photon &_old,struct photon &_out,struct crash &_count)
{         
          
          double peu_layer1 = thickness1*((ma_1+ms_1)/mut_0);
          
                         if(line_number==0)
                         {         //讀到資料第一筆資料 
                                   _old.x=_new.x;
                                   _old.y=_new.y;
                                   _old.z=_new.z;
                                   _out.x=_new.x;
                                   _out.y=_new.y;
                                   _out.z=_new.z;
                                   //fout.write((char*)&_out,sizeof _out);                         
                                   }
                         else{     //開始shift
                                  shift(peu_layer1,_old,_new,_out,_count,ma_1,ms_1,ma_2,ms_2);                                                                       
                                  _old.x=_new.x;
                                  _old.y=_new.y;
                                  _old.z=_new.z; 
                                  }
                         line_number+=1;
}

double simpson(double a, double b, int n,double Collect2Source,double r)
{
	double c=0.0,h3=0.0,result=0.0;

	

	c=(a+b)/2.0;

	h3= abs(b-a)/6.0;

	result= h3*(p(a,Collect2Source,r)+4.0*p(c,Collect2Source,r)+p(b,Collect2Source,r));

	//cout<<setw(14)<<p(a,Collect2Source,r)<<setw(14)<<p(c,Collect2Source,r)<<setw(14)<<p(b,Collect2Source,r)<<endl;

	return result;
}
double p(double x,double s,double r)
{
        double P=0.0,FirstItem=0.0,SecondItem=0.0,ThirdItem=0.0;
        
        FirstItem=(s-x);

        SecondItem=acos((s*s+(s-x)*(s-x)-illumination_r*illumination_r)/(2*(s-x)*s)) ;
		if((s*s+(s-x)*(s-x)-illumination_r*illumination_r)/(2*(s-x)*s)>=1)
		{SecondItem=0.0;}

        ThirdItem=acos((r*r+(s-x)*(s-x)-collect_r*collect_r)/(2*(s-x)*r));

		if((r*r+(s-x)*(s-x)-collect_r*collect_r)/(2*(s-x)*r)>=1)
		{ThirdItem=0.0;}
        
	    P=FirstItem*SecondItem*ThirdItem;
		//cout<< setw(14)<<ThirdItem<< setw(14)<<(r*r+(s-x)*(s-x)-collect_r*collect_r)/(2*(s-x)*r)<<endl;


       return P;
}

