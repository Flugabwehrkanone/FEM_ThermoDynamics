#pragma once

#include<cuda_runtime.h>
#include<device_launch_parameters.h>
#include<vector>
#include<Eigen/Dense>

//matrix multiplier algorithm and its kernel 
__global__ void matrixScalarMultiply_Kernel(double* matrix, double scalar, unsigned int rows, unsigned int cols);
void  matrixScalarMultiply(Eigen::MatrixXd& matrix, double scalar);


//partitioner algorithm as cuda, probably faster than the original one in utilities.h
__global__ void partitionerCuda_Kernel(double* device_vector, double start, double step, unsigned int pieces);
void partitionerCuda(std::vector<double>& vector, double start, double end, unsigned int piece);

__global__ void matrixAdd_Kernel(const double* A, const double* B, double* C, unsigned int rows, unsigned int cols);
void matrixAdd(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, Eigen::MatrixXd& C);

__global__ void matrixVectorMultiply_Kernel(const double* matrix, const double* vector, double* device_result, unsigned int rows, unsigned int cols);
Eigen::VectorXd matrixVectorMultiply(const Eigen::MatrixXd& matrix, const Eigen::VectorXd& vector);
