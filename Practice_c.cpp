// Practice_c.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include"utilities.h"
#include"Matrices.h"
#include"Solver.h"
#include<chrono>
#include<fstream>
#include"matrix_ops.cuh"

int main()
{
	auto start = std::chrono::high_resolution_clock::now();

	unsigned short int poly_degree = 5;
	unsigned int n = 2;	
	long double tau = 0.3;
	long double kappa = 0.000008L;
	long double dt = 0.001;
	long double time_start = 0.0;
	long double time_end = 10.0;
	double lambda = 3;
	
	//Eigen::MatrixXd a =alfa_solver(poly_degree, n, dt, kappa, tau);
	
	auto end = std::chrono::high_resolution_clock::now();
	/*std::ofstream alfa("alfa.txt");
	alfa << a.row(0);
	alfa.close();*/
	std::cout << E_matrix(poly_degree, n, dt, 1, 0.01, kappa, lambda, tau) << std::endl;
	
	std::chrono::duration<double> duration = end - start;
	std::cout << "Function took " << duration.count() << " seconds." << std::endl;
	//std::cout << a.block(0,0,(poly_degree * n - n + 1), ceil((time_end - time_start) / dt)) << std::endl;
	
	return 0;
	
}

