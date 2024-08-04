#include"Matrices.h"

//GENERAL UPDATE: all matrices with integration has been updated to Simpson's rule, accuracy has increased
//GENERAL UPDATE: All matrices with integration has been updated to Gaussian Quadrature, accuracy at its peak, further optimization might be required 
//All matrices are giving the correct values compared to Octave

//K matrix(N*dN with integration)

Eigen::MatrixXd K_matrix(unsigned short int poly_degree) {
	Eigen::MatrixXd K = Eigen::MatrixXd::Zero(poly_degree, poly_degree);
	// Numerical integration to populate the matrix

#pragma omp parallel for collapse(2)
	for ( short int row = 0; row < poly_degree; row++) {
		for ( short int col = 0; col < poly_degree; col++) {
			K(row, col) = GaussianQuadrature(integrated_Leg, integrated_Leg_d, row + 1, col + 1);
			}
		}
	return K;
}

//N*N matrix(N*N with integration)

Eigen::MatrixXd NN_matrix(unsigned short int poly_degree) {
	Eigen::MatrixXd NN = Eigen::MatrixXd::Zero(poly_degree, poly_degree);
	// Numerical integration to populate the matrix
#pragma omp parallel for collapse(2)
	for (short int row = 1; row < poly_degree + 1; row++) {
		for (short int col = 1; col < poly_degree + 1; col++) {
			NN(row - 1, col - 1) = GaussianQuadrature(integrated_Leg, integrated_Leg, row, col);
		}
	}
	return NN;
}

//dN*dN matrix(dN*dN with integration)

Eigen::MatrixXd dNdN_matrix(unsigned short int poly_degree) {
	Eigen::MatrixXd K = Eigen::MatrixXd::Zero(poly_degree, poly_degree);
	// Numerical integration to populate the matrix
#pragma omp parallel for collapse(2)
	for (short int row = 1; row < poly_degree + 1; row++) {
		for (short int col = 1; col < poly_degree + 1; col++) {
			K(row - 1, col - 1) = GaussianQuadrature(integrated_Leg_d, integrated_Leg_d, row, col);
		}
	}
	return K;
}

//N*dN matrix without integration, works perfectly

Eigen::MatrixXd NdN_matrix(unsigned short int poly_degree, long double length) {
	Eigen::MatrixXd M = Eigen::MatrixXd::Zero(poly_degree, poly_degree);
	for (unsigned short int row = 0; row < poly_degree; row++) {
		for (unsigned short int col = 0; col < poly_degree; col++) {
			M(row, col) = integrated_Leg(row + 1, length) * integrated_Leg_d(col + 1, length);
		}
	}
	return M;
}

//A matrix, works perfectly

Eigen::MatrixXd A_matrix(unsigned short int poly_degree, unsigned int n, long double h, double lambda, double tau) {
	//cons declaration
	const short unsigned int rho = 2600;
	const short unsigned int cv = 800;
	unsigned int size = 2 * (poly_degree * n - n + 1);
	unsigned int offset = size / 2;
	unsigned int poly_1 = poly_degree - 1;
	double C_multiplier = (h / 2) * (-1 * rho * cv);
	double T_multiplier = (h / 2) * (tau / lambda);

	//matrix declaration, with sizes
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(size, size);
	Eigen::MatrixXd C = Eigen::MatrixXd::Zero(poly_degree, poly_degree);

	//C matrix into NN matrix
	C = NN_matrix(poly_degree);

	//Row and col swap for the C matrix
	C.row(1).swap(C.row(poly_1));
	C.col(1).swap(C.col(poly_1));


	//At this point  there is no difference between the C and T matrices, so we can declare the T matrix as the copy of the C matrix
	Eigen::MatrixXd T=C;

	//Row and col swap for the T matrix in not needed, as it has already been done on the C matrix

	//multiplying both matrices
	matrixScalarMultiply(C, C_multiplier);
	matrixScalarMultiply(T, T_multiplier);
	
	/*C *= C_multiplier;
	T *= T_multiplier;*/

	
	//Filling the sparse band matrix A
	for (unsigned int i = 0; i < n; i++) {
		A.block(0+i*poly_1, 0 + i * poly_1, poly_degree, poly_degree) += C;
		A.block(offset + i * poly_1, offset + i *poly_1, poly_degree, poly_degree) += T;
	}
	return A;
}

//B matrix, works perfectly, upgrading the numerical integration may enchance accuracy

Eigen::MatrixXd B_matrix(unsigned short int poly_degree, unsigned int n, long double h, double lambda, double kappa) {

	//Declaration of the matrices
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(2 * (poly_degree * n - n + 1), 2 * (poly_degree * n - n + 1));
	Eigen::MatrixXd K = K_matrix(poly_degree);
	Eigen::MatrixXd K_t(K.transpose());
	Eigen::MatrixXd I1 = dNdN_matrix(poly_degree);
	Eigen::MatrixXd I2 = NN_matrix(poly_degree);
	Eigen::MatrixXd I0 = NdN_matrix(poly_degree, -1);
	Eigen::MatrixXd In = NdN_matrix(poly_degree, 1);
	Eigen::MatrixXd I = Eigen::MatrixXd::Zero(poly_degree, poly_degree);
	double I_10n_multiplier = (kappa / lambda) * (2.0L / h);
	double I2_multiplier = (1.0L / lambda * (h / 2.0L));
	unsigned short int poly_1 = poly_degree - 1;

	//definition of the matrices, row and col swap, multiplying the matrices with constants
	K.row(1).swap(K.row(poly_1));
	K.col(1).swap(K.col(poly_1));
	K_t.row(1).swap(K_t.row(poly_1));
	K_t.col(1).swap(K_t.col(poly_1));
	
	matrixScalarMultiply(I1, I_10n_multiplier);
	matrixScalarMultiply(I2, I2_multiplier);
	matrixAdd(I1, I2, I);
	I.row(1).swap(I.row(poly_1));
	I.col(1).swap(I.col(poly_1));

	matrixScalarMultiply(I0, I_10n_multiplier);
	I0.row(1).swap(I0.row(poly_1));
	I0.col(1).swap(I0.col(poly_1));

	matrixScalarMultiply(In, I_10n_multiplier);
	In.row(1).swap(In.row(poly_1));
	In.col(1).swap(In.col(poly_1));

	unsigned int offset = n * poly_degree - n + 1;

	//Filling the B matrix with the submatrices
	for (unsigned int i = 0; i < n; i++) {
		B.block(offset + i * poly_1, 0 + i * poly_1, poly_degree, poly_degree) += K;
		B.block(0 + i * poly_1, offset + i * poly_1, poly_degree, poly_degree) += K_t;
		B.block(offset + i * poly_1,offset + i * poly_1, poly_degree, poly_degree) += I;
	}
	B.block(offset, offset, poly_degree, poly_degree) += I0;
	B.bottomRightCorner(poly_degree, poly_degree) -= In;

	return B;
}

//D matrix works perfectly

Eigen::MatrixXd D_matrix(unsigned short int poly_degree, unsigned int n, long double dt, double theta, long double h, double kappa, double lambda, double tau) {

	//Matrix declarations
	Eigen::MatrixXd A = A_matrix(poly_degree, n, h, lambda, tau);
	Eigen::MatrixXd B = B_matrix(poly_degree, n, h, lambda, kappa);
	double B_multiplier = dt * theta;
	unsigned int size = 2 * (poly_degree * n - n + 1);
	Eigen::MatrixXd D = Eigen::MatrixXd::Zero(size, size);

	matrixScalarMultiply(B, B_multiplier);
	matrixAdd(A, B, D);
	/*B *=B_multiplier ;
	D = A + B;*/
	
	return D;
}

//E matrix works perfectly

Eigen::MatrixXd E_matrix(unsigned short int poly_degree, unsigned int n, long double dt, double theta, long double h, double kappa, double lambda, double tau) {

	//Matrix declarations
	Eigen::MatrixXd A = A_matrix(poly_degree, n, h, lambda, tau);
	Eigen::MatrixXd B = B_matrix(poly_degree, n, h, lambda, kappa);
	double B_multiplier= dt * (1 - theta);
	unsigned int size = 2 * (poly_degree * n - n + 1);
	Eigen::MatrixXd E = Eigen::MatrixXd::Zero(size, size);


	matrixScalarMultiply(B, B_multiplier);
	matrixAdd(A, B, E);
	/*B *= dt * (1-theta);
	E = A - B;*/
	return E;
}