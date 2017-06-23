#include "header.cuh"
//============================================================================================================
//Monte carlo Simulation core code
//============================================================================================================

//Set up morphology of each layer
void get_media_morphology(media *dMorphology)
{
	media *hMorphology = (media *)malloc((numberOfMedium + 2)*sizeof(media));

	hMorphology[0].Thickness = 0.0;
	hMorphology[0].n = n_source;
	hMorphology[0].Absorb_of_Baseline = 0.0;
	hMorphology[0].Scattering_of_Baseline = 0.0;
	hMorphology[0].Anisotropy = 0.9;

	hMorphology[1].Thickness = TotalThickness;
	hMorphology[1].n = n_Medium;
	hMorphology[1].Absorb_of_Baseline = 0.0;
	hMorphology[1].Scattering_of_Baseline = 200.0;
	hMorphology[1].Anisotropy = 0.9;

	hMorphology[2].Thickness = 0.0;
	hMorphology[2].n = n_Medium;
	hMorphology[2].Absorb_of_Baseline = 0.0;
	hMorphology[2].Scattering_of_Baseline = 0.0;
	hMorphology[2].Anisotropy = 0.9;

	cudaMemcpy(dMorphology, hMorphology, sizeof(media)*(numberOfMedium + 2), cudaMemcpyHostToDevice);
	free(hMorphology);
}

//Set the parameters of source fiber and detector
void get_fiber_parameters(fibers *dParameters)
{
	fibers *hParameters = (fibers *)malloc((numberOfDetectors + 1) *sizeof(fibers));
	hParameters[0].radius = illumination_r;
	hParameters[0].NA = NAOfSource;
	hParameters[0].position = 0.0;
	hParameters[0].angle = AngleOfSource * M_PI / 180;

	for (int Number = 1; Number <= numberOfDetectors; Number++) {
		hParameters[Number].radius = collect_r;
		hParameters[Number].NA = NAOfDetector;
		hParameters[Number].angle = AngleOfDetector * M_PI / 180;
		//hParameters[Number].position = ((illumination_r)+(collect_r)*(2 * Number - 1)) / cos(AngleOfDetector * M_PI / 180);
		//hParameters[Number].position = ((illumination_r)+(collect_r)*(0.1 * Number - 1)) / cos(AngleOfDetector * M_PI / 180);
		hParameters[Number].position = (0.02 + 0.02*(Number-1))/ cos(AngleOfDetector * M_PI / 180);
	}

	cudaMemcpy(dParameters, hParameters, sizeof(fibers)*(numberOfDetectors + 1), cudaMemcpyHostToDevice);
	free(hParameters);
}

__global__ void initCurand(curandState *state, unsigned long seed, int length)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx<length)	curand_init(seed, threadIdx.x, Threshold, &state[idx]);
}

//check for direction of photon whether smaller than NA of Detector or not
__device__ bool DetectionOfPhoton(fibers *FiberParameter, float ux, float uz, float Ndetector, int NumDetectors)
{
	int detected = NumDetectors + 1;
	bool detector_identified = false;
	float uz_rotated;

	for (int i = 1; i <= NumDetectors && !detector_identified; i++) {
		uz_rotated = (ux*sin(FiberParameter[i].angle)) + (uz*cos(FiberParameter[i].angle));

		if (acos(fabs(uz_rotated)) <= asin(FiberParameter[i].NA / Ndetector)) {
			detector_identified = true;
		}
	}

	return detector_identified;
}

//initial all the photon record
__global__ void initializeData(ScalingResult *data, int Lambda, int NumberOfDetectors)
{
	int gid = blockIdx.x * blockDim.x + threadIdx.x;

	if (gid < Lambda) {
		data[0].r[gid] = 0.f;
		data[0].weight[gid] = 0.f;
		data[1].r[gid] = 0.f;
		data[1].weight[gid] = 0.f;
		data[2].r[gid] = 0.f;
		data[2].weight[gid] = 0.f;
		for (int i = 0; i < NumberOfDetectors; i++) {
			data[0].r_vs_weight[gid][i] = 0.f;
			data[1].r_vs_weight[gid][i] = 0.f;
			data[2].r_vs_weight[gid][i] = 0.f;
		}
	}
}

__global__ void initializePATH(float3 *Path, int threshold)
{
	int gid = blockIdx.x * blockDim.x + threadIdx.x;

	if (gid < threshold) {
			Path[gid].x = 0.f;
			Path[gid].y = 0.f;
			Path[gid].z = 0.f;
	}
}

__global__ void testRG(float *dst, curandState *state, int length)
{
	int idx = blockDim.x*blockIdx.x + threadIdx.x;
	if (idx < length) {
		dst[idx] = curand_uniform(&state[idx]);
	}
}
//Initialize position and direction of light with consideration of NA and angle of fiber
__device__ void Launch(curandState &state, photon &P, fibers *dF, long id, float Nsrc, int Detect) //初始化光子一開始的位置、方向、權重 
{
	float rnd_position  = curand_uniform(&state);
	float rnd_Azimuth   = curand_uniform(&state);
	float rnd_direction = curand_uniform(&state);
	float rnd_rotated   = curand_uniform(&state);

	float AzimuthAngle  = 2 * M_PI*rnd_Azimuth;
	float rotated_angle = 2 * M_PI*rnd_rotated;

	float launchPosition = dF[0].radius*sqrtf(rnd_position);
	P.x = launchPosition*cos(AzimuthAngle) / cos(dF[0].angle);
	P.y = launchPosition*sin(AzimuthAngle);
	P.z = 0.0;

	float theta_direction = asin(dF[0].NA / Nsrc)*rnd_direction;

	P.ux = (sin(theta_direction)*cos(rotated_angle)*cos(dF[0].angle))
		- (cos(theta_direction)*sin(dF[0].angle));

	P.uy = sin(theta_direction)*sin(rotated_angle);

	P.uz = (sin(theta_direction)*cos(rotated_angle)*sin(dF[0].angle))
		+ (cos(theta_direction)*cos(dF[0].angle));
	if (Detect == 0) {
		P.x = 0, P.y = 0, P.z = 0;
	}
	P.weight = 1.0f;
}

__device__ float fresnel(float n_in, float n_out, float uz)
{
	float theta_in = acos(fabs(uz));
	float theta_t = asin((n_in / n_out)*sin(theta_in));

	float R = 0.5 * ((sin(theta_in - theta_t)*sin(theta_in - theta_t))
		/ (sin(theta_in + theta_t)*sin(theta_in + theta_t)) + (tan(theta_in - theta_t)*tan(theta_in - theta_t))
		/ (tan(theta_in + theta_t)*tan(theta_in + theta_t)));

	return R;
}

__device__ bool Specular(curandState &state, float uz, media *dM)
{
	bool specular_reflection = false;	// Default value is 0; If internal reflection will	occur, becomes 1 
	float theta_c;						// Critical angle for incidence 
	float rnd;							// Fresnel reflectance 
	float r;
	float n_in = 0, n_out = 0;

	n_in = dM[0].n;
	n_out = dM[1].n;

	// Total internal reflection possible
	theta_c = (n_out < n_in) ? asin(n_out / n_in) : M_PI / 2.0;


	if (acos(fabs(uz)) >= theta_c) {	// total internal reflection       
		specular_reflection = true;
	} else if (acos(fabs(uz))<theta_c) {
		r = fresnel(n_in, n_out, uz);
		rnd = curand_uniform(&state);	// Generates a random number in (0, 1) 
		if (rnd <= r) {					// Internal reflection
			specular_reflection = true;
		}
	}
	return specular_reflection;
}

__device__ void move(float3* dPATH, float *Grid, photon &P, float stepLength, int tid, int NUM)
{
	P.x = P.x + stepLength*P.ux;
	P.y = P.y + stepLength*P.uy;
	P.z = P.z + stepLength*P.uz;

	RecordPath(dPATH, Grid, P, tid, NUM);
}

//__ballot() and __clz()
__device__ int sign(float input)
{
	if (input>0)
		return 1;
	else if (input<0)
		return -1;
	else
		return 0;
}

__device__ void changeDirectionDueToRefraction(float &ux, float &uy, float &uz, media *dM)
{
	float n_in = 0, n_out = 0;

	if (uz>0) {
		n_in = dM[0].n;
		n_out = dM[1].n;
	} else if (uz <= 0) {
		n_in = dM[1].n;
		n_out = dM[0].n;
	}

	// Incident and transmitted angles, respectively
	float theta_i = acos(fabs(uz));
	float theta_t = asin((n_in*sin(theta_i)) / n_out);

	ux = (ux)*(n_in / n_out);
	uy = (uy)*(n_in / n_out);
	uz = (sign(uz))*cos(theta_t);
}

__device__ void ModifyPositionAndDirection(float3 *dPATH, float *Grid, curandState devState, photon &P, int &layer, media *dM, bool TIF, int tid, int NUM)
{
	float DIB = 0.f;

	if (P.uz<0)	{
		DIB = -P.z / P.uz;
	} else if (P.uz>0) {
		DIB = (dM[1].Thickness - P.z) / P.uz;
	}

	if (TIF) {
//		photon_record.Array_site -= 1;
		move(dPATH, Grid, P, getUnitStep(devState, layer, dM), tid, NUM);
		DIB = -DIB;
		P.uz = -P.uz;
		move(dPATH, Grid, P, getUnitStep(devState, layer, dM), tid, NUM);
	} else {
//		photon_record.Array_site -= 1;
		move(dPATH, Grid, P, getUnitStep(devState, layer, dM), tid, NUM);
		changeDirectionDueToRefraction(P.ux, P.uy, P.uz, dM);
	}

}

__device__ bool checkForTIF(curandState &state, photon P, int &layer, media *dM)
{
	int photonOut = 0;	// Default value is 0; If internal reflection will	occur, becomes 1
	float theta_c;		// Critical angle for incidence 
	float rnd;			// Fresnel reflectance 

	float r;
	float n_in = 0, n_out = 0;
	bool TIF = false;
	if (P.uz>0)	{
		n_in = dM[1].n;
		n_out = dM[2].n;
	} else if (P.uz<0) {
		n_in = dM[1].n;
		n_out = dM[0].n;
	}
	// Calculates the critical angle 
	if (n_out<n_in)	{// Total internal reflection possible
		theta_c = asin(n_out / n_in);

		
		if (acos(fabs(P.uz)) >= theta_c) {	// Total internal reflection
			TIF = true;
		} else	{
			r = fresnel(n_in, n_out, P.uz);
			rnd = curand_uniform(&state);	// Generates a random number in (0, 1) 
			// Internal reflection
			if (rnd <= r)	TIF = true;
		}
	} else {
		r = fresnel(n_in, n_out, P.uz);
		rnd = curand_uniform(&state);	// Generates a random number in (0, 1) 
		// Internal reflection
		if (rnd <= r)	TIF = true;
	}
	return TIF;
}



__device__ float InitializationOfWeight(int layer, media *dM, float Nsrc)
{
	float R_sp = ((Nsrc - dM[layer].n)*(Nsrc - dM[layer].n)) / ((Nsrc + dM[layer].n)*(Nsrc + dM[layer].n));
	return 1 - R_sp;	//return the weight
}

__device__ float getUnitStep(curandState &state, int layer, media *dM)
{
	float rn = 0.0;
	do{
		rn = curand_uniform(&state);
	} while (rn == 0);

	return (-1 * log(rn) / (float)(dM[1].Absorb_of_Baseline + dM[1].Scattering_of_Baseline));
}


__device__ void direction(curandState &state, float &ux, float &uy, float &uz, int layer, media *dM)
{
	float rn, temp, theta, phi;
	float uxprime, uyprime, uzprime;
	float g = dM[1].Anisotropy;

	if (g == 0) {
		rn = curand_uniform(&state);
		theta = acos(2 * rn - 1);
	} else {
		rn = curand_uniform(&state);
		temp = (1 - g*g) / (1 - g + 2 * g*rn);
		theta = acos((1 / (2 * g))*(1 + g*g - (temp*temp)));
	}

	rn = curand_uniform(&state);

	phi = 2 * M_PI * rn;

	if (abs(uz)>0.9999) {
		uxprime = sin(theta)*cos(phi);
		uyprime = sin(theta)*sin(phi);
		uzprime = uz*cos(theta) / abs(uz);
	} else {
		uxprime = sin(theta) / sqrt(1 - uz*uz)*(ux*uz*cos(phi) - uy*sin(phi)) + ux*cos(theta);
		uyprime = sin(theta) / sqrt(1 - uz*uz)*(uy*uz*cos(phi) + ux*sin(phi)) + uy*cos(theta);
		uzprime = -sin(theta)*cos(phi)*sqrt(1 - ((uz)*(uz))) + uz*cos(theta);
	}

	ux = uxprime;
	uy = uyprime;
	uz = uzprime;
}

//判斷 光子是否存活  (當光子被吸收到很小時才會需要做此判斷 (在test_in_tissue會用到)) 
__device__ void roulette(bool &absorb, float &weight, unsigned int seed, photon P, int layer, media *dM, float roulette_m)
{
	curandState state;
	curand_init(seed, 0, 0, &state);
	float rn = curand_uniform(&state);
	absorb = false;
	if (rn<(1 / roulette_m)) {
		weight = roulette_m*weight;
		direction(seed, P.ux, P.uy, P.uz, layer, dM);
	} else {
		weight = 0;
		absorb = true;
	}
}

__device__ void absorption(float &WA, float &weight, int layer, media *dM)
{
	WA = WA + weight*(dM[layer].Absorb_of_Baseline / (dM[layer].Absorb_of_Baseline + dM[layer].Scattering_of_Baseline));
	weight = weight*(dM[layer].Scattering_of_Baseline / (dM[layer].Absorb_of_Baseline + dM[layer].Scattering_of_Baseline));
}

__device__ void newPosition(float &z, float &uz, int layer, media *dM)
{
	if (z<0) {
		z = -z;
		uz = -uz;
	} else if (z >dM[layer].Thickness) {
		z = 2 * dM[layer].Thickness - z;
		uz = -uz;
	}
}


__device__ void RecordPath(float3 *dPATH, float *Grid, photon P, int tid, int NUM)
{
	if (NUM < Threshold) {
		//recorded path for scaling
		dPATH[NUM + tid * Threshold].x = P.x;
		dPATH[NUM + tid * Threshold].y = P.y;
		dPATH[NUM + tid * Threshold].z = P.z;

		//record path for building the distribution of photons
		int i, j, k;
		int LR = sizeR / deltaR;
		int LZ = sizeZ / deltaZ;
		unsigned int idx;
		//if (NUM > 0)
		{
			i = int(P.x / deltaR) + (LR / 2);
			j = int(P.y / deltaR) + (LR / 2);
			k = int(P.z / deltaZ);

			if (i < LR && j < LR && k < LZ && i >= 0 && j >= 0 && k >= 0)
			{
				idx = (int)(i + j * LR + k * LR * LR);
				Grid[idx] +=  (1 - (NUM / Threshold));
			}
		}
	}
}

__device__ void ClearRecordPath(float3 *dPATH, int tid)
{
	unsigned int idx = 0;
	for (int i = 0; i < Threshold; i++)
	{
		idx = tid * Threshold + i;
		dPATH[idx].x = 0.f;
		dPATH[idx].y = 0.f;
		dPATH[idx].z = 0.f;
	}
}

//Load optical parameter of two layer from "input_parameter_two_layer.txt" file
//=========================================================================================================================================================================================================================================
void get_optical_property(float & UpLayerCoefficient, float & A_UpLayer, float & K_UpLayer, float & Hb_BottomLayer, float & Sto2_BottomLayer, float & A_BottomLayer, float &K_BottomLayer, float &Thickness_UpLayer, float &Collagen)
{
	ifstream finOpticalProperty;
	finOpticalProperty.open("input_parameter_two_layer.txt", ios::in);

	if (!finOpticalProperty) { cout << "fail read file" << endl; }

	finOpticalProperty >> UpLayerCoefficient >> A_UpLayer >> K_UpLayer >> Hb_BottomLayer >> Sto2_BottomLayer >> A_BottomLayer >> K_BottomLayer >> Thickness_UpLayer >> Collagen;
	A_UpLayer = A_UpLayer * 10000;
	A_BottomLayer = A_BottomLayer * 10000;
}

//輸入 extinct coefficient 並計算 光學參數
void mua_data(float *mua_out, float f, float sto2, int length, float col)
{
	float  lamda, hemo, oxhemo, collagen_ab;
	int i = 0;
	int j = 0;

	ifstream infile;
	infile.open("epsilon.txt");
	
	ifstream col_para;
	col_para.open("collagen_mua_data.txt");

	while (infile >> lamda >> oxhemo >> hemo) {
		mua_out[i] = f*((sto2)*2.303*oxhemo / 64532 + (1 - sto2)*2.303*hemo / 64500);
		i++;
	}

	while (col_para >> lamda >> collagen_ab){
		mua_out[j] = mua_out[j] + col*collagen_ab;
		j++;
	}
	infile.close();
	col_para.close();
}

void mua_data_up(float *mua_out, float up_c, int length)
{
	float a = 0.0, b = 0.0;
	int i = 0;

	FILE *fin = fopen("epsilon_2.txt", "r");

	while (fscanf(fin, "%f %f", &a, &b) != EOF)	{
		mua_out[i] = up_c*b;
		i++;
	}
	fclose(fin);
}

void mus_data(float *mus, float A, float K, int length)
{
	int j = 0;

	for (float i = 250; i <= 1000; i += 2) {
		mus[j] = A*pow(i, -K);
		j++;
	}
}
//=============================================================================================================
void DataOutput(ScalingResult *data)
{
	ofstream foutForTotalRefelectance;
	ofstream foutForReflectanceWithDistance;

	foutForTotalRefelectance.open("C://Users//AZSXDCFVGBHNJ4Y6//Desktop//SMC_GPU//x64//Release//two_layer_mc_weight_with_scale.txt");
	foutForReflectanceWithDistance.open("C://Users//AZSXDCFVGBHNJ4Y6//Desktop//SMC_GPU//x64//Release//two_layer_mc_weight_with_scale_and_r.txt");

	int x = 0;

	for (int l = 250; l <= 1000; l += 2) {
		foutForTotalRefelectance << l << "\t" << data[0].weight[x] / NumberOfPhotons << "\n";
		foutForReflectanceWithDistance << l << "\t";

		for (int i = 0; i<numberOfDetectors; i++) {
			foutForReflectanceWithDistance << fixed << data[0].r_vs_weight[x][i] / NumberOfPhotons << "\t";
		}

		foutForReflectanceWithDistance << "\n";
		x++;
	}
};


__global__ void runPhoton(photon *P, ScalingResult *dData, curandState *devState, float3 *dPATH, float *Grid, float *ma1, float *ma2, float *ms1, float *ms2,
	float Thickness_UpLayer, media *dM, fibers *dF, float Nsrc, float Ndetector, int Detect, int numPhoton, int NumDetectors)
{
	int gid = blockIdx.x * blockDim.x + threadIdx.x;
	int tid = threadIdx.x;

	if (gid < numPhoton)
	{
		Launch(devState[gid], P[gid], dF, gid, Nsrc, Detect);

		if (!Specular(devState[gid], P[gid].uz, dM))
		{
			int layer = 0;
			bool TIF = false, Absorb = false, Reflect = false, Transmit = false;
			int run_number = 0;
			changeDirectionDueToRefraction(P[gid].ux, P[gid].uy, P[gid].uz, dM);
			while (!Absorb && !Reflect && !Transmit)
			{
				layer = 1;
				move(dPATH, Grid, P[gid], getUnitStep(devState[gid], layer, dM), tid, run_number);

				if (P[gid].z < 0) {
					if (checkForTIF(devState[gid], P[gid], layer, dM)) {
						TIF = true;
						run_number++;
						ModifyPositionAndDirection(dPATH, Grid, devState[gid], P[gid], layer, dM, TIF, tid, run_number);
						direction(devState[gid], P[gid].ux, P[gid].uy, P[gid].uz, layer, dM);						
					} else {
						TIF = false;
						ModifyPositionAndDirection(dPATH, Grid, devState[gid], P[gid], layer, dM, TIF, tid, run_number);
						Reflect = true;
					}
				} else if (P[gid].z > dM[1].Thickness) {					
					if (checkForTIF(devState[gid], P[gid], layer, dM)) {
						TIF = true;
						run_number++;
						ModifyPositionAndDirection(dPATH, Grid, devState[gid], P[gid], layer, dM, TIF, tid, run_number);
						direction(devState[gid], P[gid].ux, P[gid].uy, P[gid].uz, layer, dM);						
					} else {
						TIF = false;
						Transmit = true;
					}
				} else {
					direction(devState[gid], P[gid].ux, P[gid].uy, P[gid].uz, layer, dM);
					run_number++;
				}

				if (run_number > Threshold) { ClearRecordPath(dPATH, tid); Absorb = true; } //視為吸收

				if (Reflect && run_number < Threshold &&  DetectionOfPhoton(dF, P[gid].ux, P[gid].uz, Ndetector, NumDetectors)) {
					Scaling(dData, dPATH, tid, run_number, ma1, ms1, ma2, ms2, Thickness_UpLayer, dM, dF);
				}
			}			
			ClearRecordPath(dPATH, tid);
			//P[gid].z = dPATH[blockSize * Threshold - 1].z;// mut_0(dM) / (ms1[275] + ma1[275]);
		}		
	}
	__syncthreads();
}

//==========================================================================================================
__device__ float mut_0(media *dM)
{
	return dM[1].Absorb_of_Baseline + dM[1].Scattering_of_Baseline;
}

__device__ float afa_0(media *dM)
{
	return dM[1].Scattering_of_Baseline / (dM[1].Absorb_of_Baseline + dM[1].Scattering_of_Baseline);
}

__device__ float mut_1(float ms_1, float ma_1)	{ return ms_1 + ma_1 ; }
__device__ float afa_1(float ms_1, float mut_1)	{ return ms_1 / mut_1; }
__device__ float mut_2(float ms_2, float ma_2)	{ return ms_2 + ma_2 ; }
__device__ float afa_2(float ms_2, float mut_2)	{ return ms_2 / mut_2; }


__device__ void Scaling(ScalingResult *dData, float3 *dPATH, int tid, int NUM, float *ma1_out, float *ms1_out, float *ma2_out, float *ms2_out, float Thickness_UpLayer, media *dM, fibers *dF)
{
	for (int i = 75; i<276; i++)
	{		
		int2 count_crash;
		count_crash.x = 0;
		count_crash.y = 0;
		int line_number = 0;
		int CounterOfPhotons = 0;
		int run_number = NUM;

		float3 newP;
		float3 oldP;
		float3 outP;
		newP.x = 0, newP.y = 0, newP.z = 0;
		oldP.x = 0, oldP.y = 0, oldP.z = 0;
		outP.x = 0, outP.y = 0, outP.z = 0;

		float ma_1 = ma1_out[i], ms_1 = ms1_out[i], ma_2 = ma2_out[i], ms_2 = ms2_out[i];
		float mut0 = mut_0(dM), mut1 = ms_1 + ma_1, mut2 = ms_2 + ma_2;
		float c1 = mut0 / mut1, c2 = mut0 / mut2;
		float peu_layer = Thickness_UpLayer / c1;
		bool Flag = true;
		
		do{
			if (run_number == CounterOfPhotons)
			{
				Flag = false;
			}
			else
			{
				newP = dPATH[CounterOfPhotons + tid * Threshold];
				scale(c1, c2, peu_layer, line_number, newP, oldP, outP, count_crash);
				CounterOfPhotons++;
			}
		} while (Flag && run_number < CounterOfPhotons);
		//dPATH[blockSize * Threshold - 1].z = outP.x;

		float r = (outP.x * outP.x + outP.y * outP.y > 0.0f) ? sqrtf(outP.x * outP.x + outP.y * outP.y) : 0.f;

		dData[0].r[i] = r;

		float weight = 1 * (pow(afa_1(ms_1, mut1) / afa_0(dM), count_crash.x - 1))*(pow(afa_2(ms_2, mut2) / afa_0(dM), count_crash.y));

		dData[0].weight[i] += weight;

		float FirstItem = 0.0, SecondItem = 0.0, ThirdItem = 0.0;

		if (detect == 0)
		{
			float Collect2Source = 0.0;
			float lo = 0.0, up = 0.0;
			float p[10] = { 0.f };

			for (int j = 0; j<numberOfDetectors; j++)
			{
				p[j] = 0.0;
				Collect2Source = dF[j + 1].position;
				lo = (-illumination_r > Collect2Source - r - collect_r) ? -illumination_r : Collect2Source - r - collect_r;
				up = (illumination_r  < Collect2Source - r + collect_r) ? illumination_r  : Collect2Source - r + collect_r;

				if (Collect2Source >= (illumination_r + collect_r) && (r >= (Collect2Source - illumination_r - collect_r)) && (r<(Collect2Source + illumination_r + collect_r)))
				{
					for (float x = lo; x<up; x += 0.00005)
					{
						float FirstItem = 0.0, SecondItem = 0.0, ThirdItem = 0.0;
						FirstItem = (Collect2Source - x);

						SecondItem = acos((Collect2Source*Collect2Source + (Collect2Source - x)*(Collect2Source - x) - illumination_r*illumination_r) / (2 * (Collect2Source - x)*Collect2Source));
						if ((Collect2Source*Collect2Source + (Collect2Source - x)*(Collect2Source - x) - illumination_r*illumination_r) / (2 * (Collect2Source - x)*Collect2Source) >= 1)
						{
							SecondItem = 0.0;
						}

						ThirdItem = acos((r*r + (Collect2Source - x)*(Collect2Source - x) - collect_r*collect_r) / (2 * (Collect2Source - x)*r));

						if ((r*r + (Collect2Source - x)*(Collect2Source - x) - collect_r*collect_r) / (2 * (Collect2Source - x)*r) >= 1)
						{
							ThirdItem = 0.0;
						}

						p[j] += FirstItem*SecondItem*ThirdItem;
					}
					p[j] *= (2 / (M_PI*M_PI*illumination_r*illumination_r))*0.00005;
				}
				else
				{
					p[j] = 0;
				}
				if (p[j] <= FLT_MAX && p[j] >= -FLT_MAX)
				{
					dData[0].r_vs_weight[i][j] += (weight*p[j]);
				}
				else
				{
					dData[0].r_vs_weight[i][j] += 0;
				}
			}
		}
		else if (detect == 1)
		{
			for (int f = 1; f <= numberOfDetectors; f++)
			{
				if ((((outP.x - dF[f].position) / (dF[f].radius / cos(dF[f].angle)))*
					((outP.x - dF[f].position) / (dF[f].radius / cos(dF[f].angle)))) +
					(((outP.y) / dF[f].radius)*((outP.y) / dF[f].radius)) <= 1.0)
				{
					dData[0].r_vs_weight[i][f - 1] += weight;
				}
			}
		}
	}
}

__device__ void scale(float c1, float c2, float peu_layer1, int &line_number, float3 &newP, float3 &oldP, float3& outP, int2 &crash)
{
	//float peu_layer1 = (c1 == 0.f) ? thickness1 : thickness1 / c1;
	c1 = 1, c2 = 1;
	if (line_number == 0) {         //讀到資料第一筆資料 
		oldP = newP;
		outP = newP;
	} else{     //開始shift
		
		float3 p;

		if (oldP.z <= peu_layer1)
		{
			crash.x++;
			if (newP.z <= peu_layer1)
			{
				outP.x += c1*(newP.x - oldP.x);
				outP.y += c1*(newP.y - oldP.y);
				outP.z += c1*(newP.z - oldP.z);
			}
			else	//===================================================================
			{                //定義在peusudo layer的位置
				p.x = (newP.z - oldP.z == 0.f) ? oldP.x : oldP.x + (peu_layer1 - oldP.z) / (newP.z - oldP.z)*(newP.x - oldP.x);
				p.y = (newP.z - oldP.z == 0.f) ? oldP.x : oldP.y + (peu_layer1 - oldP.z) / (newP.z - oldP.z)*(newP.y - oldP.y);
				p.z = peu_layer1;
				//===================================================================
				//step3
				//縮放在每一層pseudo layer內的水平位移，加總起來得到總共的水平位移
				outP.x += (c1*(p.x - oldP.x) + c2*(newP.x - p.x));
				outP.y += (c1*(p.y - oldP.y) + c2*(newP.y - p.y));
				outP.z += (c1*(p.z - oldP.z) + c2*(newP.z - p.z));
			}
		}
		else
		{
			crash.y++;
			if (newP.z >= peu_layer1)
			{
				outP.x += c2*(newP.x - oldP.x);
				outP.y += c2*(newP.y - oldP.y);
				outP.z += c2*(newP.z - oldP.z);
			}
			else
			{                //定義在peusudo layer的位置
				p.x = (newP.z - oldP.z == 0.f) ? oldP.x : oldP.x + (peu_layer1 - oldP.z) / (newP.z - oldP.z)*(newP.x - oldP.x);
				p.y = (newP.z - oldP.z == 0.f) ? oldP.y : oldP.y + (peu_layer1 - oldP.z) / (newP.z - oldP.z)*(newP.y - oldP.y);
				p.z = peu_layer1;
				//=================================================================
				//step3
				//縮放在每一層pseudo layer內的水平位移，加總起來得到總共的水平位移
				outP.x += (c1*(newP.x - p.x)+c2*(p.x - oldP.x));
				outP.y += (c1*(newP.y - p.y)+c2*(p.y - oldP.y));
				outP.z += (c1*(newP.z - p.z)+c2*(p.z - oldP.z));
				//==================================================================
			}
		}
		oldP = newP;
	}
	line_number++;
}

__device__ float simpson(float a, float b, int n, float Collect2Source, float r)
{
	float c = (a + b) / 2.0;
	float h3 = abs(b - a) / 6.0;
	float result = h3*(p(a, Collect2Source, r) + 4.0*p(c, Collect2Source, r) + p(b, Collect2Source, r));

	return result;
}

__device__ float p(float x, float s, float r)
{
	float FirstItem = 0.0, SecondItem = 0.0, ThirdItem = 0.0;

	FirstItem = (s - x);

	float temp1 = (s*s + (s - x)*(s - x) - illumination_r*illumination_r) / (2 * (s - x)*s);
	SecondItem = (temp1 >= 1) ? 0.f : acos(temp1);

	float temp2 = (r*r + (s - x)*(s - x) - collect_r*collect_r) / (2 * (s - x)*r);
	ThirdItem = (temp2 >= 1) ? 0.f : acos(temp2);

	return FirstItem*SecondItem*ThirdItem;
}

//==========================================================================================================
void MonteCarlo(photon * dPhotonArray, ScalingResult *dData, float *Grid, float *ma1, float *ma2, float *ms1, float *ms2
	, float &thickness_layer, media *dM, fibers *dF, float Nsrc, float Ndetector, int Detect, int PhotonNumber)
{
	//define the size of block and grid
	int gridSize = (PhotonNumber + blockSize - 1) / blockSize;
	int NumDetector = numberOfDetectors;
	unsigned long seed = time(NULL);
	
	// 'PATH' for recording the path of photons (1024 * Steps = blockSize * Threshold)
	float3 *hPATH = (float3 *)malloc(blockSize*Threshold*sizeof(float3));
	float3 *dPATH;
	cudaMalloc((void**)&dPATH, sizeof(float3)* blockSize*Threshold);
	initializePATH << <Threshold, blockSize >> >(dPATH, blockSize*Threshold);

	// initial the state of random number on CUDA
	curandState *devState;
	cudaMalloc((void**)&devState, PhotonNumber * sizeof(curandState));
	initCurand << <gridSize,blockSize >> >(devState, seed, PhotonNumber);
	
	// initial the Scaling Results
	int Lambda = LamdaNumber;
	initializeData << <ceil(Lambda / blockSize), blockSize >> >(dData, Lambda, NumDetector);
	
	// Start the Monte Carlo with scaling
	runPhoton << <gridSize, blockSize >> >(dPhotonArray, dData, devState, dPATH, Grid, ma1, ma2, ms1, ms2,
		thickness_layer, dM, dF, Nsrc, Ndetector, Detect, PhotonNumber, NumDetector);

	cudaMemcpy(hPATH, dPATH, sizeof(float3)*blockSize*Threshold, cudaMemcpyDeviceToHost);
	cudaFree(devState);
	cudaFree(dPATH);


	// validate the coordinates of photon on the 'PATH' array
	ofstream PATHFILE;
	PATHFILE.open("C://Users//AZSXDCFVGBHNJ4Y6//Desktop//SMC_GPU//x64//Release//PATH.txt");

	for (int i = 0; i < Threshold; i++) {
		PATHFILE << i << ":\t" << hPATH[i].x << "\t" << hPATH[i].y << "\t" << hPATH[i].z << "\t" << hPATH[1023 * Threshold + i].x << "\t" << hPATH[1023 * Threshold + i].y << "\t" << hPATH[1023 * Threshold + i].z << endl;
	}
	PATHFILE.close();
	free(hPATH);
}


void MC_Migraiton()
{
	// info. of media on device
	media *dM;
	cudaMalloc((void **)&dM, sizeof(media)*(numberOfMedium + 2));
	get_media_morphology(dM);

	// info. of fiber on device
	fibers *dF;
	cudaMalloc((void **)&dF, sizeof(fibers)*(numberOfDetectors + 1));
	get_fiber_parameters(dF);

	// setting the coefficients of absorbing and scattering
	float UpLayerCoefficient, A_UpLayer, K_UpLayer, Hb_BottomLayer, Sto2_BottomLayer, A_BottomLayer, K_BottomLayer, Thickness_UpLayer, Collagen;
	get_optical_property(UpLayerCoefficient, A_UpLayer, K_UpLayer, Hb_BottomLayer, Sto2_BottomLayer, A_BottomLayer, K_BottomLayer, Thickness_UpLayer, Collagen);

	float *h_ma1_out = (float *)malloc(LamdaNumber*sizeof(float));
	float *h_ma2_out = (float *)malloc(LamdaNumber*sizeof(float));
	float *h_ms1_out = (float *)malloc(LamdaNumber*sizeof(float));
	float *h_ms2_out = (float *)malloc(LamdaNumber*sizeof(float));
	mua_data_up(h_ma1_out, UpLayerCoefficient               , LamdaNumber);		//get上層吸收係數 
	mus_data   (h_ms1_out, A_UpLayer      , K_UpLayer       , LamdaNumber);		//get上層散射係數 
	mua_data   (h_ma2_out, Hb_BottomLayer , Sto2_BottomLayer, LamdaNumber, Collagen);		//get下層吸收係數
	mus_data   (h_ms2_out, A_BottomLayer  , K_BottomLayer   , LamdaNumber);		//get下層散射係數 

	float *d_ma1_out, *d_ma2_out, *d_ms1_out, *d_ms2_out;
	cudaMalloc((void **)&d_ma1_out, sizeof(float)*LamdaNumber);
	cudaMalloc((void **)&d_ma2_out, sizeof(float)*LamdaNumber);
	cudaMalloc((void **)&d_ms1_out, sizeof(float)*LamdaNumber);
	cudaMalloc((void **)&d_ms2_out, sizeof(float)*LamdaNumber);
	cudaMemcpy(d_ma1_out, h_ma1_out, sizeof(float)*LamdaNumber, cudaMemcpyHostToDevice);
	cudaMemcpy(d_ma2_out, h_ma2_out, sizeof(float)*LamdaNumber, cudaMemcpyHostToDevice);
	cudaMemcpy(d_ms1_out, h_ms1_out, sizeof(float)*LamdaNumber, cudaMemcpyHostToDevice);
	cudaMemcpy(d_ms2_out, h_ms2_out, sizeof(float)*LamdaNumber, cudaMemcpyHostToDevice);
	free(h_ma1_out);
	free(h_ms1_out);
	free(h_ma2_out);
	free(h_ms2_out);

	// allocate the memory to Scaling Results
	ScalingResult *hData = (ScalingResult *)malloc(3 * sizeof(ScalingResult));
	ScalingResult *dData;
	cudaMalloc((void **)&dData , sizeof(ScalingResult)*3);

	// allocate the memory to grid structure
	int GridR = sizeR / deltaR;
	int GridZ = sizeZ / deltaZ;
	float *hGrid = (float *)malloc(GridR * GridR * GridZ * sizeof(float));
	float *dGrid;
	cudaMalloc((void **)&dGrid, sizeof(float)* GridR * GridR * GridZ);
	cudaMemset(dGrid, 0.0f, GridR * GridR * GridZ * sizeof(float));

	// allocating the memory to record all the photon info.	
	int defaultNumPhoton = ceil((float)NumberOfPhotons / 1024.0) * 1024;
	cout << "Photon Number     : " << defaultNumPhoton << endl;
	cout << "<<<Memory Requirement>>>" << endl;
	cout << "Photon Array      : " << (float)defaultNumPhoton / 1024 / 1024 * sizeof(photon) << "(MB)" << endl;	
	cout << "Photon's PATH     : " << (float)blockSize*Threshold / 1024 / 1024 * sizeof(float3) << "(MB)" << endl;
	cout << "Optical Parameters: " << (float)LamdaNumber * 4 / 1024 / 1024 * sizeof(float) << "(MB)" << endl;
	cout << "Scaling Result    : " << (float)1 / 1024 / 1024 * sizeof(ScalingResult) << "(MB)" << endl;
	cout << "Grid Structure    : " << (float)GridR * GridR * GridZ / 1024 / 1024 * sizeof(float) << "(MB)" << endl;

	photon *hPhotoArray = (photon *)malloc(defaultNumPhoton * sizeof(photon));
	photon *dPhotonArray;
	cudaMalloc((void**)&dPhotonArray, sizeof(photon)* defaultNumPhoton);
		

	float Nsrc = n_source;
	float Ndetector = n_detector;
	int Detect = detect;

	MonteCarlo(dPhotonArray, dData, dGrid, d_ma1_out, d_ma2_out, d_ms1_out, d_ms2_out, Thickness_UpLayer, dM, dF, Nsrc, Ndetector, Detect, defaultNumPhoton);
	cudaMemcpy(hPhotoArray, dPhotonArray, sizeof(photon)*defaultNumPhoton, cudaMemcpyDeviceToHost);
	cudaMemcpy(hGrid, dGrid, sizeof(float)*GridR * GridR * GridZ, cudaMemcpyDeviceToHost);
	cudaMemcpy(hData, dData, sizeof(ScalingResult)* 3, cudaMemcpyDeviceToHost);

	// export the Scaling Results (hData) and Grid Structure
	DataOutput(hData);

	FILE *fp;
	fp= fopen("C://Users//AZSXDCFVGBHNJ4Y6//Desktop//SMC_GPU//x64//Release//test.raw", "wb+");
	if (fp)	{
		fwrite(hGrid, GridR * GridR * GridZ, sizeof(float), fp);
		fclose(fp);
	} else { printf("Export Failed."); }
		
	// review the final coordinates of latest photon (should not be 0,0,0) 
	printf("<<TEST>> P%d: (%f, %f, %f) \n", defaultNumPhoton - 1, hPhotoArray[defaultNumPhoton - 1].x, hPhotoArray[defaultNumPhoton - 1].y, hPhotoArray[defaultNumPhoton-1].z);


	cudaFree(dM);
	cudaFree(dF);
	cudaFree(d_ma1_out);
	cudaFree(d_ms1_out);
	cudaFree(d_ma2_out);
	cudaFree(d_ms2_out);
	cudaFree(dData);
	cudaFree(dPhotonArray);
	cudaFree(dGrid);
	free(hData);
	free(hPhotoArray);
	free(hGrid);
}
//==================================================================