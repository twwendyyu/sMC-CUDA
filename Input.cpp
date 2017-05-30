#include "header.h"
using namespace std;

//Load optical parameter of two layer from "input_parameter_two_layer.txt" file


//=========================================================================================================================================================================================================================================
void get_optical_property(float & UpLayerCoefficient,float & A_UpLayer,float & K_UpLayer,float & Hb_BottomLayer,float & Sto2_BottomLayer,float & A_BottomLayer,float &K_BottomLayer,float &Thickness_UpLayer)
{
     
    ifstream finOpticalProperty;
    
    finOpticalProperty.open("input_parameter_two_layer.txt",ios::in);

	if(!finOpticalProperty)
	{
		cout<<"fail read file"<<endl;

	}

    
    finOpticalProperty>>UpLayerCoefficient>>A_UpLayer>>K_UpLayer>>Hb_BottomLayer>>Sto2_BottomLayer>>A_BottomLayer>>K_BottomLayer>>Thickness_UpLayer;

	A_UpLayer=A_UpLayer*10000;

	A_BottomLayer=A_BottomLayer*10000;
    
    
    return;
    
}

//輸入 extinct coefficient 並計算 光學參數
struct mua_array mua_data(float f,float sto2)
{
	float  lamda [LamdaNumber];
    float  hemo  [LamdaNumber];
    float  oxhemo[LamdaNumber];
    mua_array mua_out;
	for(int i=0;i<LamdaNumber;i++)
	{
		lamda[i]=0.0;
		hemo[i]=0.0;   
		oxhemo[i]=0.0;    
		mua_out.all[i]=0.0;
	}

	int i=0;

	ifstream infile;
    infile.open("epsilon.txt");
	
	while(infile>>lamda[i]>>oxhemo[i]>>hemo[i])
	 {
		   
           mua_out.all[i]=f*((sto2)*2.303*oxhemo[i]/64532+(1-sto2)*2.303*hemo[i]/64500);
           i++;
	 }
	infile.close();
	
	return mua_out;
 }

 struct mua_array mua_data_up(double up_c)
{
    float lamda[LamdaNumber];
    float mua_up[LamdaNumber];
    float a=0.0,b=0.0,c=0.0,d=0.0;  
	mua_array mua_out_up;
	for(int i=0;i<LamdaNumber;i++)
	{
		lamda[i]=0.0;
		mua_up[i]=0.0;
		mua_out_up.all[i]=0.0;
	}
  
	int i=0;

	FILE *fin = fopen("epsilon_2.txt", "r");

	 while(fscanf(fin,"%f %f", &a, &b)!= EOF)
	 {
           lamda[i]=a;
           mua_up[i]=b;
           mua_out_up.all[i]=up_c*mua_up[i];
           i++;
	 }
	
	fclose(fin);
	
	return mua_out_up;
 }

struct mus_array mus_data(float A,float K)
{       
	 
    int j=0,a=0;

	mus_array mus;

	for (float i=250;i<=1000;i+=2)
	{
		mus.all[j]=A*pow(i,-K);
		j++;
	}

	return mus;
 }
//=============================================================================================================
    
    
    

//=========================================================================================================================================================================================================================================

//Load extinction coefficent from "epsilon.txt" file

//=========================================================================================================================================================================================================================================
/*std::vector<vector<double > > get_extinction_Coeffiecint()
{
     std::vector<vector<double > >ExtincCoef(3);

	 ifstream finExtinctionCoefficent;
	 
     finExtinctionCoefficent.open("epsilon.txt",ios::in);
     
     if(!finExtinctionCoefficent)
     {
         cout<<"epsilon file read fail"<<endl;

		 exit(-1);

     }
     std::string temp;
     
     getline(finExtinctionCoefficent, temp, '\x1a');
     
     std::istringstream sinput(temp);
     
     double  i;
     
	 while(sinput>>i)
	 {           
           ExtincCoef[0].push_back(i);
           
           sinput>>i;
           
           ExtincCoef[1].push_back(i);
           
           sinput>>i;
           
           ExtincCoef[2].push_back(i);       
	 }

	return ExtincCoef;
 }*/
//=========================================================================================================================================================================================================================================
