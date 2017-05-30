#include "header.h"

using namespace std;

// parameters of fibers ,morphology of media

struct ScalingResult data;
struct ScalingResult data2;
struct ScalingResult data3;

extern struct fibers parameters[numberOfDetectors+1];

extern struct media morphology[numberOfMedium+2];

extern struct ArrayOfPhotonPath photon_record;

struct mua_array mua_data(float f,float sto2);
struct mus_array mus_data(float A,float K);
struct mua_array mua_data_up(double up_c);

//===========================================================================================================

void RecordPath(double &,double & ,double &);

//All the function of photon propagration
void Launch(double &x,double &y,double &z,double &ux,double &uy,double &uz);
void checkExitLayer(double &z,double &uz,int &layer);
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
void Scaling(struct  mua_array &ma1_out,struct  mus_array &ms1_out,struct  mua_array &ma2_out,struct  mus_array &ms2_out,float &Thickness_UpLayer);
bool DetectionOfPhoton(double &ux, double &uz);
//============================================================================================================

//============================================================================================================
void  get_fiber_parameters();
void  get_media_morphology();
void get_optical_property(float & ,float & ,float & ,float & ,float & ,float & ,float &,float &);
void initializeData();

double random(long *idum);
//============================================================================================================

//Monte carlo Simulation core code

//============================================================================================================

long seed;
//long int tim;
__time64_t tim;

void initializeData()
{ 
	photon_record.Array_site=0;

	 for(int j=0;j<LamdaNumber;j++)
	 {
		 data.weight[j]=0.0;
		 data.r[j]=0.0;
		 for(int n=0;n<numberOfDetectors;n++)
		 {
			 data.r_vs_weight[j][n]=0;
			 data2.r_vs_weight[j][n]=0;
			 data3.r_vs_weight[j][n]=0;
		 }
	 } 

	 for(int s1=0;s1<3;s1++)
	 {
		 for(int s=0;s<Threshold;s++)
		 {
			 photon_record.Array[s1][s]=0;
		 }
	 }

}

int MC_Migraiton()
{
	      srand((unsigned int)time(NULL));
		  /*seed=tim=0;
		  seed=abs((long) time(&tim));
		  
		  long dummy;
		  dummy=-1;
		  random(&dummy);*/

	      double x=0,y=0,z=0;
	      double ux=0,uy=0,uz=0;
		  double WA=0,WR=0,WT=0;
		  double weight=0;
		  double unitStep=0;
		  bool   TIF=false;
		  int count=0;

		  float UpLayerCoefficient, A_UpLayer,K_UpLayer, Hb_BottomLayer,Sto2_BottomLayer,A_BottomLayer,K_BottomLayer,Thickness_UpLayer;
		 
		  bool Absorb,Transmit,Reflect;
		  int photonOut;
		  int run_number,layer;

		  get_fiber_parameters();
		  get_media_morphology();

		  get_optical_property(UpLayerCoefficient, A_UpLayer,K_UpLayer, Hb_BottomLayer,Sto2_BottomLayer,A_BottomLayer,K_BottomLayer,Thickness_UpLayer);

		  mua_array ma1_out;
		  
		  ma1_out =  mua_data_up(UpLayerCoefficient); //get糷Μ玒计 
		  
		  mus_array ms1_out;
		  
		  ms1_out = mus_data(A_UpLayer,K_UpLayer);    //get糷床甮玒计 
		  
		  mua_array ma2_out;
		  
		  ma2_out = mua_data(Hb_BottomLayer,Sto2_BottomLayer); //get糷Μ玒计
		  
		  mus_array ms2_out;
		  
		  ms2_out = mus_data(A_BottomLayer,K_BottomLayer);    //get糷床甮玒计 

		  initializeData();
		 /* ofstream aa;
		  aa.open("aa.txt", ios::app);
		  ofstream bb;
		  bb.open("bb.txt", ios::app);*/
	          
		 

		  for(int i=0;i<NumberOfPhotons;i++)
		  {
			  Launch(x,y,z,ux,uy,uz); 

              layer=0;

			  RecordPath(x,y,z);

              if(!Specular(uz))
			  {	
				  weight=1;       

				  Absorb=false,Reflect=false,Transmit=false,photonOut=0,run_number=0;	

				  changeDirectionDueToRefraction(ux,uy,uz,layer);

				  while(!Absorb && !Reflect && !Transmit)
				  {
					  layer=1;
					  getUnitStep(unitStep,layer);	  
					  move(x,y,z,ux,uy,uz,unitStep);  
					  if (z<0)
					  {
						  if(checkForTIF(x,y,z,ux,uy,uz,layer,weight,WR,WA,WT))
						  {
							  TIF=true;
							  ModifyPositionAndDirection(x,y,z,ux,uy,uz,layer,TIF);		  
							  direction(ux,uy,uz,layer); 
							  run_number++;
						  }
						  else
						  {
							  TIF=false;
							  ModifyPositionAndDirection(x,y,z,ux,uy,uz,layer,TIF);
							  Reflect=true;
						  }
					  }
					  else if (z>morphology[1].Thickness)
					  {   
						  if(checkForTIF(x,y,z,ux,uy,uz,layer,weight,WR,WA,WT))
						  {
							  TIF=true;
							  ModifyPositionAndDirection(x,y,z,ux,uy,uz,layer,TIF);
							  direction(ux,uy,uz,layer);  
							  run_number++;
						  }
						  else
						  {
							  TIF=false;
							  Transmit=true;
						  }
					  }
					  else
					  {
						  direction(ux,uy,uz,layer);
						  run_number++;
					  }					  
					  
					  if(run_number>Threshold){ClearRecordPath();Absorb=true;} //跌Μ
					  if(Reflect && run_number<=Threshold &&  DetectionOfPhoton(ux,uz))
					  {
						  count++;
						  Scaling(ma1_out,ms1_out,ma2_out,ms2_out,Thickness_UpLayer);
						  //cout<<count<<"\n";
						  //bb<<photon_record.Array_site<<endl;
						 /* for(int q=0;q<photon_record.Array_site;q++)
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
				   
				  ClearRecordPath();	
              } 		  	   
		  }	    
		  cout<<count<<"\n";
		  //cout<<WR/NumberOfPhotons<<endl;
		  return 0;
}
