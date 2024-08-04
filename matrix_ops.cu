#include"matrix_ops.cuh"



//////////////////////////////////////////////////		Multiplying matrix with scalar		/////////////////////////////////////////////////////////////////////////////////

__global__ void matrixScalarMultiply_Kernel(double* matrix, double scalar, unsigned int rows, unsigned int cols) {
	int row = blockIdx.x * blockDim.x + threadIdx.x;
	int col = blockIdx.y * blockDim.y + threadIdx.y;

	if (row < rows && col < cols) {
		int index = row * cols + col;
		matrix[index] *= scalar;
	}
}

void matrixScalarMultiply(Eigen::MatrixXd& matrix, double scalar) {
	int rows = matrix.rows();
	int cols = matrix.cols();
	size_t size = rows * cols * sizeof(double);

	double* device_matrix;
	cudaMalloc(&device_matrix, size);
	cudaMemcpy(device_matrix, matrix.data(), size, cudaMemcpyHostToDevice);

	dim3 threadsPerBlock(16, 16);
	dim3 numBlocks((rows + threadsPerBlock.x - 1) / threadsPerBlock.x, (cols + threadsPerBlock.y - 1) / threadsPerBlock.y);

	matrixScalarMultiply_Kernel <<<numBlocks, threadsPerBlock >>> (device_matrix, scalar, rows, cols); //intellisense gives error for the triple <<< and >>>,ignore it, code can run

	cudaDeviceSynchronize();

	cudaMemcpy(matrix.data(), device_matrix, size, cudaMemcpyDeviceToHost);

	cudaFree(device_matrix);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////		The partitioner algorithm		/////////////////////////////////////////////////////////////////////////////////

__global__ void partitionerCuda_Kernel(double* device_vector, double start, double step, unsigned int pieces) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < pieces) {
		device_vector[idx] = start + idx * step;
	}
}

void partitionerCuda(std::vector<double>& vector, double start, double end, unsigned int piece) {
	vector.resize(piece);
	if (start > end) {
		std::swap(start, end);
	}
	double step = (end - start) / (piece - 1);

	double* device_vector;
	size_t size = piece * sizeof(double);
	cudaMalloc(&device_vector, size);
	unsigned short int threadsPerBlock = 256;
	unsigned int numBlocks = (piece + threadsPerBlock - 1) / threadsPerBlock;

	partitionerCuda_Kernel <<<numBlocks, threadsPerBlock >>> (device_vector, start, step, piece);

	cudaDeviceSynchronize();
	cudaMemcpy(vector.data(), device_vector, size, cudaMemcpyDeviceToHost);
	cudaFree(device_vector);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////	 Matrix addition kernel and function		/////////////////////////////////////////////////////////////

__global__ void matrixAdd_Kernel(const double* A, const double* B, double* C, unsigned int rows, unsigned int cols) {
	unsigned int row = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int col = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int index = row * cols + col;
	if (row < rows && col < cols) {
		C[index] = A[index] + B[index];
	}
}
void matrixAdd(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, Eigen::MatrixXd& C) {
	
	unsigned int rows = A.rows();
	unsigned int cols = A.cols();
	size_t size = rows * cols * sizeof(double);

	double* device_A;
	double* device_B;
	double* device_C;

	cudaMalloc(&device_A, size);
	cudaMalloc(&device_B, size);
	cudaMalloc(&device_C, size);

	cudaMemcpy(device_A, A.data(), size, cudaMemcpyHostToDevice);
	cudaMemcpy(device_B, B.data(), size, cudaMemcpyHostToDevice);

	dim3 threadsPerBlock(16, 16);
	dim3 numBlocks((rows + threadsPerBlock.x - 1) / threadsPerBlock.x, (cols + threadsPerBlock.y - 1) / threadsPerBlock.y);

	matrixAdd_Kernel <<<numBlocks, threadsPerBlock >>> (device_A, device_B, device_C, rows, cols);

	cudaDeviceSynchronize();
	cudaMemcpy(C.data(), device_C, size, cudaMemcpyDeviceToHost);

	cudaFree(device_A);
	cudaFree(device_B);
	cudaFree(device_C);

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////		Multiplying matrix with vector		/////////////////////////////////////////////////////////////////////////

__global__ void matrixVectorMultiply_Kernel(const double* matrix, const double* vector, double* device_result, unsigned int rows, unsigned int cols) {

	int row = blockIdx.x * blockDim.x + threadIdx.x;

	if (row < rows) {
		double sum = 0.0;
		for (unsigned int col = 0; col < cols; col++) {
			sum += matrix[row * cols + col] * vector[col];
		}
		device_result[row] = sum;
	}
}

Eigen::VectorXd matrixVectorMultiply(const Eigen::MatrixXd& matrix,const Eigen::VectorXd& vector) {

	unsigned int rows = matrix.rows();
	unsigned int cols = matrix.cols();
	Eigen::VectorXd result(rows);

	size_t sizeMatrix = rows * cols * sizeof(double);
	size_t sizeVector = rows * sizeof(double);

	double* device_matrix;
	double* device_vector;
	double* device_result;

	cudaMalloc(&device_matrix, sizeMatrix);
	cudaMalloc(&device_vector, sizeVector);
	cudaMalloc(&device_result, sizeVector);

	cudaMemcpy(device_matrix, matrix.data(), sizeMatrix, cudaMemcpyHostToDevice);
	cudaMemcpy(device_vector, vector.data(), sizeVector, cudaMemcpyHostToDevice);

	int blockSize = 256;
	int gridSize = (rows + blockSize - 1) / blockSize;

	matrixVectorMultiply_Kernel <<< gridSize, blockSize >>> (device_matrix, device_vector, device_result, rows, cols);

	cudaMemcpy(result.data(), device_result, sizeVector, cudaMemcpyDeviceToHost);

	cudaFree(device_matrix);
	cudaFree(device_vector);
	cudaFree(device_result);

	return result;

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
