#pragma once
#include<vector>
#include<cmath>
#include<stdexcept>
#include<iostream>
#include<functional>
#include<Eigen/Dense>

//Polynomial functions, all tested, all working fine
 long double Legendre(unsigned short int poly_degree, long double  x);
 long double Legendre_d(unsigned short int poly_degree, long double x);
// void Legendre_d_helper(unsigned short int poly_degree, long double x, long double& Pn, long double& Pn_1);
 long double integrated_Leg(unsigned short int poly_degree, long double x);
 long double integrated_Leg_d(unsigned short int poly_degree, long double x);

 std::vector<long double>Function_integr_values(std::function<long double(unsigned short int, long double)> func1, std::function<long double(unsigned short int, long double)> func2,
	unsigned int poly_degree1, unsigned int poly_degree2);

 long double Simpsons_rule(const std::vector<long double>& functionValues, long double a, long double b);

//Partitioner function, creates interval with given start and end points interval lenght is piece+1
std::vector<long double> partitioner(long double interval_start, long double interval_end, unsigned int piece);

//Col and row swapping functions for std::vector matrices, not needed after replacing std::vector funciotns with eigen matrices
//void swapColumns(std::vector<std::vector<long double>>& matrix, unsigned short int col1, unsigned short int col2);
//void Matrix_multiply(std::vector<std::vector<long double>>& matrix, float multiplier);


//Test functions, working fine, not implemented yet to other functions, or to the main, may need to rewrite them to eigen vectors
//UPDATE: they are rewritten to eigen::vectors
Eigen::VectorXd Test_function(long double time, unsigned short int poly_degree, unsigned int n);
Eigen::VectorXd Laser_heating(long double time, unsigned short int poly_degree, unsigned int n);
Eigen::VectorXd f_test(unsigned short int poly_degree, unsigned int n, double theta, long double time, long double dt);

