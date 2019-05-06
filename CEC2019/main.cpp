/*
  CEC19 Test Function Suite for Single Objective Optimization
  Noor Awad (email: noorawad1989@gmail.com) 
  Sep. 21th 2018
*/

#include <iomanip>
#include <stdio.h>
#include <malloc.h>
#include "cec19_func.hpp"

int main (int argc, char *argv[]) {
	Cec19 cec19 = Cec19();

	double *f, *x;
	char FileName[48];
	int i, j, k, n, m, func_num;

	FILE *fpt;
	m = 2;
	n = 10;
	x = (double *)malloc(m * n *sizeof(double));
	f = (double *)malloc(sizeof(double)  *  m);

	for (i = 0; i < 10; i++) {
		func_num = i + 1;
		sprintf(FileName, "input_data/shift_data_%d.txt", func_num);

		fpt = fopen(FileName, "r");
		if (fpt == NULL) {
			std::cout << "Error: Cannot open input file for reading" << std::endl;
		}

		if (x == NULL) {
			std::cout << "Error: there is insufficient memory available" << std::endl;
		}

		for (k = 0; k < n; k++) {
			fscanf(fpt, "%lf", &x[k]);
		}

		fclose(fpt);

		for (j = 0; j < n; j++) {
			x[1 * n + j] = 0.0;
		}

		for (k = 0; k < 1; k++) {
			cec19.test_func(x, f, n, m, func_num);
			for (j = 0; j < 2; j++) {
				std::cout << "-> f" << func_num << "(x[" << (j + 1) << "]) = " << std::setprecision(16) << f[j] << std::endl;
			}

			std::cout << std::endl;
		}
	}

	return 1;
}
