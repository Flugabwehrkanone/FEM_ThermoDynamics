#pragma once

#include<vector>
#include"utilities.h"
#include<omp.h>
#include<Eigen/Dense>
#include<Eigen/Core>

//small matrices (p*p)
//UPDATE: all matrices are rewritten from std::vector to eigen matrices

Eigen::MatrixXd K_matrix(unsigned short int poly_degree); 
Eigen::MatrixXd NN_matrix(unsigned short int poly_degree);
Eigen::MatrixXd dNdN_matrix(unsigned short int poly_degree);
Eigen::MatrixXd NdN_matrix(unsigned short int poly_degree, long double length);

//Large matrices (n*p-n+1)*2
Eigen::MatrixXd A_matrix(unsigned short int poly_degree, unsigned int n, long double h, double lambda, double tau);
Eigen::MatrixXd B_matrix(unsigned short int poly_degree, unsigned int n, long double h, double lambda, double kappa);
Eigen::MatrixXd D_matrix(unsigned short int poly_degree, unsigned int n, long double dt, double theta, long double h, double kappa, double lambda, double tau);
Eigen::MatrixXd E_matrix(unsigned short int poly_degree, unsigned int n, long double dt, double theta, long double h, double kappa, double lambda, double tau);
