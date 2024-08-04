#pragma once

#include"Matrices.h"
#include"matrix_ops.cuh"

Eigen::MatrixXd alfa_solver(unsigned short int poly_degree, unsigned int n, long double dt, double kappa, double tau);