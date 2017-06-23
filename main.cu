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
	system("pause");



	return 0;
}