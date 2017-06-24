#pragma once
#include "header.cuh"


int main()
{
	clock_t start;   //Timer
	long double diff;
	start = clock();

	MC_Migraiton();

	diff = (std::clock() - start) / (long double)CLOCKS_PER_SEC;
	cout << "computation time: " << diff << '\n';
	
	
	#if defined(_WIN32) || defined(_WIN64)
		system("pause");
	#endif

	return 0;
}