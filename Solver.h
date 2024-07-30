#pragma once
#include<Eigen/Dense>
#include"Matrices.h"
#include<functional>
#include<omp.h>

Eigen::MatrixXd alfa_solver(unsigned short int poly_degree, unsigned int n, long double dt, double kappa, double tau);