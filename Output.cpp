#include "header.h"
using namespace std;


extern struct ScalingResult data;
extern struct ScalingResult data2;
extern struct ScalingResult data3;
extern struct fibers parameters[numberOfDetectors+1];

void get_fiber_parameters();

//Data Output

void DataOutput()
{
	ofstream foutForTotalRefelectance;
	ofstream foutForReflectanceWithDistance;

	foutForTotalRefelectance.open("two_layer_mc_weight_with_scale.txt");
	foutForReflectanceWithDistance.open("two_layer_mc_weight_with_scale_and_r.txt");
	/*fin2.open("top.txt");
	fin3.open("bottom.txt");*/
	get_fiber_parameters();

	int x=0;

	 for (int l=250;l<=1000;l+=2)          
	 {                               
		 foutForTotalRefelectance << l <<"\t"<< data.weight[x]/NumberOfPhotons<<"\n";
                
		 foutForReflectanceWithDistance<< l <<"\t";
		/* fin2<<l<<"\t";
		 fin3<<l<<"\t";*/
                
		 for(int i=0;i<numberOfDetectors;i++)
                 
		 {
			 foutForReflectanceWithDistance<<fixed<<data.r_vs_weight[x][i]/NumberOfPhotons<<"\t";
			 /*fin2<<setprecision(8)<<(data2.r_vs_weight[x][i]/data.r_vs_weight[x][i])/NumberOfPhotons<<"\t";
			 fin3<<setprecision(8)<<(data3.r_vs_weight[x][i]/data.r_vs_weight[x][i])/NumberOfPhotons<<"\t";*/
		 }

		 foutForReflectanceWithDistance<<"\n";
		/* fin2<<"\n";
		 fin3<<"\n";*/
                 
		 x++;                      
	 }
};
