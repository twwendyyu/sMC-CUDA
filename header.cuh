#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <cfloat>
#include <vector>
#include <string>
#include <sstream> 
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <device_launch_parameters.h>

#define blockSize 1024

#define NumberOfPhotons 10000000
#define TotalThickness 1.0
#define Threshold  2048


//在此區調整光纖參數
#define n_source 1.61
#define n_detector 1.61
#define n_Medium 1.4
#define NAOfSource 0.4
#define NAOfDetector 0.4
//#define AngleOfSource -46.0	//toward x axis :positive
#define AngleOfSource 0 //Wang modified
//#define AngleOfDetector  46.0		//toward x axis :negative
#define AngleOfDetector 0  //Wang modified
#define detect  1						// convolution:0 another:1


#define m 20
#define w_th 0.005
#define numberOfDetectors 10
#define numberOfMedium 2
#define LamdaNumber 376
#define illumination_r 0.01
#define collect_r 0.01


#define sizeR 1.5
#define sizeZ 1.5
#define deltaR 0.005
#define deltaZ 0.005


using namespace std;

//程式應用結構
struct ArrayOfPhotonPath
{
	int Array_site;
	float Array[3][Threshold];
	float vector[2];
};

struct GridSystemForReflectance
{
	float Array[100][100];
};

struct GridSystemForAbsorbance
{
	float Array[100][100];
};

struct ScalingResult
{
	float weight[LamdaNumber];
	float r[LamdaNumber];
	float r_vs_weight[LamdaNumber][10];
};

struct crash
{
	float crash1, crash2;
};

struct photon
{
	float x, y, z, ux, uy, uz, weight;
};

struct media
{
	float Thickness;
	float n;
	float Anisotropy;
	float Absorb_of_Baseline;
	float Scattering_of_Baseline;
};

struct fibers
{
	float radius;			/* Radius of the fiber in cm */
	float NA;				/* NA of the fiber (in air) */
	float position;			/* Distance from the origin to the center position of the
							fiber in cm (defined to be 0.0 for the source fiber) */
	float angle;			/* Tilt angle of the fiber in radians */
};


void get_media_morphology(media *dMorphology);
void get_fiber_parameters(fibers *dParameters);
__global__ void initCurand(curandState *state, unsigned long seed, int length);
__device__ bool DetectionOfPhoton(fibers *FiberParameter, float ux, float uz, float Ndetector, int NumDetectors);
__global__ void initializeData(ScalingResult *data, int Lambda, int NumberOfDetectors);
__global__ void initializePATH(float3 *Path, int threshold);
__device__ void Launch(curandState &devState, photon &P, fibers *dF, long id, float Nsrc, int Detect);
__device__ float fresnel(float n_in, float n_out, float uz);
__device__ bool Specular(curandState &devState, float uz, media *dM);
__device__ void move(float3* dPATH, float *Grid, photon &P, float stepLength, int NUM);
__device__ int sign(float input);
__device__ void changeDirectionDueToRefraction(float &ux, float &uy, float &uz, media *dM);
__device__ void ModifyPositionAndDirection(float3 *dPATH, float *Grid, curandState devState, photon &P, int &layer, media *dM, bool TIF, int tid, int NUM);
__device__ bool checkForTIF(curandState &devState, photon P, int &layer, media *dM);
__device__ float InitializationOfWeight(int layer, media *dM, float Nsrc);
__device__ float getUnitStep(curandState &devState, int layer, media *dM);
__device__ float reflect(float WR, float weight);
__device__ float transmit(float WT, float weight);
__device__ void direction(unsigned int seed, float &ux, float &uy, float &uz, int layer, media *dM);
__device__ void roulette(bool &absorb, float &weight, unsigned int seed, photon P, int layer, media *dM, float roulette_m);
__device__ void absorption(float &WA, float &weight, int layer, media *dM);
__device__ void newPosition(float &z, float &uz, int layer, media *dM);
__device__ void RecordPath(float3 *dPATH, float *Grid, photon P, int tid, int NUM);
__device__ void ClearRecordPath(float3 *dPATH, int tid);
__device__ void ClearRecordPath();
void get_optical_property(float & UpLayerCoefficient, float & A_UpLayer, float & K_UpLayer, float & Hb_BottomLayer, float & Sto2_BottomLayer, float & A_BottomLayer, float &K_BottomLayer, float &Thickness_UpLayer);
void mua_data(float *mua_out, float f, float sto2, int length, float col);
void mua_data_up(float *mua_out, float up_c, int length);
void mus_data(float *mus, float A, float K, int length);
void DataOutput(ScalingResult data);
__global__ void runPhoton(photon *P, ScalingResult *dData, curandState *devState, float3 *dPATH, float *Grid, float *ma1, float *ma2, float *ms1, float *ms2,
	float Thickness_UpLayer, media *dM, fibers *dF, float Nsrc, float Ndetector, int Detect, int numPhoton, int NumDetectors);
__device__ float mut_0(media *dM);
__device__ float afa_0(media *dM);
__device__ float mut_1(float ms_1, float ma_1);
__device__ float afa_1(float ms_1, float mut_1);
__device__ float mut_2(float ms_2, float ma_2);
__device__ float afa_2(float ms_2, float mut_2);
__device__ void Scaling(ScalingResult *data, float3 *dPATH, int tid, int NUM, float *ma1_out, float *ms1_out, float *ma2_out, float *ms2_out, float Thickness_UpLayer, media *dM, fibers *dF);
__device__ void scale(float c1, float c2, float thickness1, int &line_number, float3 &newP, float3 &oldP, float3 &outP, int2 &count);
__device__ float simpson(float a, float b, int n, float Collect2Source, float r);
__device__ float p(float x, float s, float r);

void MonteCarlo(photon * dPhotonArray, ScalingResult *dData, float *Grid, float *ma1, float *ma2, float *ms1, float *ms2
	, float &thickness_layer, media *dM, fibers *dF, float Nsrc, float Ndetector, int Detect, int PhotonNumber);
void MC_Migraiton();
