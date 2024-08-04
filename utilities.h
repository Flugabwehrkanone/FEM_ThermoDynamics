#pragma once

#include<vector>
#include<cmath>
#include<stdexcept>
#include<iostream>
#include<functional>
#include<Eigen/Dense>
#include<Eigen/Core>

//Polynomial functions, all tested, all working fine
 long double Legendre(unsigned short int poly_degree, long double  x);
 long double Legendre_d(unsigned short int poly_degree, long double x);
 long double integrated_Leg(unsigned short int poly_degree, long double x);
 long double integrated_Leg_d(unsigned short int poly_degree, long double x);

//Partitioner function, creates interval with given start and end points interval lenght is piece+1
std::vector<long double> partitioner(long double interval_start, long double interval_end, unsigned int piece);

//Test functions, working fine, not implemented yet to other functions, or to the main, may need to rewrite them to eigen vectors
//UPDATE: they are rewritten to eigen::vectors
Eigen::VectorXd Test_function(long double time, unsigned short int poly_degree, unsigned int n);
Eigen::VectorXd Laser_heating(long double time, unsigned short int poly_degree, unsigned int n);
Eigen::VectorXd f_test(unsigned short int poly_degree, unsigned int n, double theta, long double time, long double dt);

//pure chatgpt code, needs to be testes, corrected, implemented, optimized...
//UPDATE:Code rewritten, and upgraded, waiting for implementation
Eigen::MatrixXd gauss_legendre(unsigned short int poly_degree);

//The main part of the Gaussian quadrature
long double GaussianQuadrature(std::function<long double(unsigned short int, long double)>func1, std::function<long double(unsigned short int, long double)>func2,
	unsigned short int poly_degree1, unsigned short int poly_degree2);