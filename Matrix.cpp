#include"Matrix.h"
#include"utilities.h"
#include<vector>
#include<iostream>

std::vector<std::vector<long double>>K_matrix(unsigned int poly_degree) {
	std::vector<std::vector<long double>>K;
	long double temp = 0.0;
	std::vector<double>integr_interval;
	partitioner(-1, 1, 10).swap(integr_interval);
	for (int row = 0; row != poly_degree; row++) {
		for (int col = 0; col != poly_degree; col++) {
			for (int i = 1; i < 11; i++) {
				temp = temp + Legendre(row + 1, integr_interval[i-1]) * Legendre_d(col + 1, integr_interval[i-1]);
			}
			K[row][col] = temp;
		}
	}
	return K;
}