#include"Matrices.h"


//GENERAL UPDATE: all matrices with integration has been updated to Simpson's rule, accuracy has increased

//K matrix(N*dN with integration)

Eigen::MatrixXd K_matrix(unsigned short int poly_degree) {
	Eigen::MatrixXd K = Eigen::MatrixXd::Zero(poly_degree, poly_degree);
	// Numerical integration to populate the matrix
	#pragma omp parallel for collapse(2)
	for ( short int row = 0; row < poly_degree; row++) {
		for ( short int col = 0; col < poly_degree; col++) {
			std::vector<long double> integr_interval = Function_integr_values(integrated_Leg, integrated_Leg_d, row+1, col+1 );
			K(row, col) = Simpsons_rule(integr_interval, -1, 1);
			}
		}
	return K;
}

//N*N matrix
//Still needs update for increasing accuracy

Eigen::MatrixXd NN_matrix(unsigned short int poly_degree) {
	Eigen::MatrixXd K = Eigen::MatrixXd::Zero(poly_degree, poly_degree);
	// Numerical integration to populate the matrix
	#pragma omp parallel for collapse(2)
	for ( short int row = 0; row < poly_degree; row++) {
		for ( short int col = 0; col < poly_degree; col++) {
			std::vector<long double> integr_interval = Function_integr_values(integrated_Leg, integrated_Leg, row + 1, col + 1);
			K(row, col) = Simpsons_rule(integr_interval, -1, 1);
		}
	}
	return K;
}

//dN*dN matrix; tested separately only works up to poly_degree=3,, dunno why yet; implemented into B_matrix, works perfectly
//update, the problem was the legendre derivative. The commented formula could not take 1 or -1 as x, because it resulted in division by 0.
//Still needs update for increasing accuracy

Eigen::MatrixXd dNdN_matrix(unsigned short int poly_degree) {
	Eigen::MatrixXd K = Eigen::MatrixXd::Zero(poly_degree, poly_degree);
	
	// Numerical integration to populate the matrix
#pragma omp parallel for collapse(2)
	for ( short int row = 0; row < poly_degree; row++) {
		for ( short int col = 0; col < poly_degree; col++) {
			std::vector<long double> integr_interval = Function_integr_values(integrated_Leg_d, integrated_Leg_d, row + 1, col + 1);
			K(row, col) = Simpsons_rule(integr_interval, -1, 1);
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

	//matrix declaration, with sizes
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(size, size);
	Eigen::MatrixXd C = Eigen::MatrixXd::Zero(poly_degree, poly_degree);
	
	//C matrix into NN matrix
	C = NN_matrix(poly_degree);


	//Row and col swap for the C matrix
	C.row(1).swap(C.row(poly_1));
	C.col(1).swap(C.col(poly_1));
	
	//At this point  there is no difference between the C and T matrices, so we can declare the T matrix as the copy of the C matrix
	Eigen::MatrixXd T(C);

	//multiplying both matrices
	C *= (h / 2) * (-1 * rho * cv);
	T *= (h / 2) * (tau / lambda);

	
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

	//definition of the matrices, row and col swap, multiplying the matrices with constants
	K.row(1).swap(K.row(poly_degree - 1));
	K.col(1).swap(K.col(poly_degree - 1));

	
	I1 *= (kappa / lambda) * (2.0L / h);
	I1.row(1).swap(I1.row(poly_degree - 1));
	I1.col(1).swap(I1.col(poly_degree - 1));

	I2 *= (1.0L / lambda * (h / 2.0L));
	I2.row(1).swap(I2.row(poly_degree - 1));
	I2.col(1).swap(I2.col(poly_degree - 1));

	I0 *= (kappa / lambda * (2.0 / h));
	I0.row(1).swap(I0.row(poly_degree - 1));
	I0.col(1).swap(I0.col(poly_degree - 1));

	In *= (kappa / lambda * (2.0 / h));
	In.row(1).swap(In.row(poly_degree - 1));
	In.col(1).swap(In.col(poly_degree - 1));

	Eigen::MatrixXd I = I1 + I2;

	unsigned int offset = n * poly_degree - n + 1;
	unsigned int blockSize = poly_degree - 1;

	//Filling the B matrix with the submatrices
	for (unsigned int i = 0; i < n; i++) {
		B.block(offset + i * blockSize, 0 + i * blockSize, poly_degree, poly_degree) += K;
		B.block(0 + i * blockSize, offset + i * blockSize, poly_degree, poly_degree) += K_t;
		B.block(offset + i * blockSize,offset + i * blockSize, poly_degree, poly_degree) += I;
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
	unsigned int size = 2 * (poly_degree * n - n + 1);
	Eigen::MatrixXd D = Eigen::MatrixXd::Zero(size, size);

	B *= dt * theta;
	D = A + B;
	return D;
}

//E matrix works perfectly

Eigen::MatrixXd E_matrix(unsigned short int poly_degree, unsigned int n, long double dt, double theta, long double h, double kappa, double lambda, double tau) {

	//Matrix declarations
	Eigen::MatrixXd A = A_matrix(poly_degree, n, h, lambda, tau);
	Eigen::MatrixXd B = B_matrix(poly_degree, n, h, lambda, kappa);
	unsigned int size = 2 * (poly_degree * n - n + 1);
	Eigen::MatrixXd E = Eigen::MatrixXd::Zero(size, size);

	B *= dt * (1-theta);
	E = A - B;
	return E;
}