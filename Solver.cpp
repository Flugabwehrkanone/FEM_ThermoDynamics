#include"Solver.h"

Eigen::MatrixXd alfa_solver(unsigned short int poly_degree, unsigned int n, long double dt, double kappa, double tau) {
	
	Eigen::setNbThreads(omp_get_max_threads()-2); // Set to maximum available threads
	std::cout << "Using " << Eigen::nbThreads() << " threads." << std::endl;

	const int size = 2 * (poly_degree * n - n + 1);
	long double length = 0.005L;
	const long double h = length / n;
	const double theta = 1.0;
	const double T0 = 293.0;
	//const double q0 = 0.0;
	const double lambda = 3.0;
	const long double time_start = 0.0L;
	const long double time_end = 10.0L;
	unsigned int time_steps = ceil((abs(time_end - time_start)) / dt);

	//Trying to use dynamic memory allocation for optimizing the code, still under development

	/*Eigen::MatrixXd* K = new Eigen::MatrixXd(poly_degree, poly_degree);
	Eigen::MatrixXd* NN = new Eigen::MatrixXd(poly_degree, poly_degree);
	Eigen::MatrixXd* dNdN = new Eigen::MatrixXd(poly_degree, poly_degree);
	Eigen::MatrixXd* NdN_I0 = new Eigen::MatrixXd(poly_degree, poly_degree);
	Eigen::MatrixXd* NdN_In = new Eigen::MatrixXd(poly_degree, poly_degree);

	*K = K_matrix(poly_degree);
	*NN = NN_matrix(poly_degree);
	*dNdN = dNdN_matrix(poly_degree);
	*NdN_I0 = NdN_matrix(poly_degree, -1);
	*NdN_In = NdN_matrix(poly_degree, 1);


	Eigen::MatrixXd* A = new Eigen::MatrixXd(size, size);
	Eigen::MatrixXd* B = new Eigen::MatrixXd(size, size);*/

	//A = A_matrix(poly_degree, n, h, lambda, tau, *NN);


	Eigen::VectorXd f = Eigen::VectorXd::Zero(size);
	Eigen::VectorXd alfa_col = Eigen::VectorXd::Zero(size);
	Eigen::MatrixXd E = E_matrix(poly_degree, n, dt, theta, h, kappa, lambda, tau);
	Eigen::MatrixXd D = D_matrix(poly_degree, n, dt, theta, h, kappa, lambda, tau);
	Eigen::MatrixXd alfa = Eigen::MatrixXd::Zero(size, time_steps);

	Eigen::FullPivLU<Eigen::MatrixXd> luD(D);
	Eigen::MatrixXd Dinv = luD.inverse();
	//Before the for loop something slows the code down, have to find out
	//UPDATE:the inverse of the D matirx has slowed the code down, it is a little faster now

	//alfa(0, 0) = T0;
#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		alfa(i+i*(poly_degree-1), 0) = T0;
	}
	//std::cout << alfa.col(0);

	//#pragma omp parallel for private(alfa_col, f)
	//Should make a vector matrix cuda multiplication function, to speed this for loop up
	//UPDATE: the cuda multiplication for vector and matrix is implemented
	for (int i = 1; i < time_steps; i++) {
		alfa_col = alfa.col(i - 1);
		f = f_test(poly_degree, n, theta, time_start + (i - 1) * dt, dt);
		alfa.col(i) = (f + matrixVectorMultiply(E, alfa_col));
		alfa.col(i) = matrixVectorMultiply(Dinv, alfa.col(i));
		//std::cout << i << std::endl;
	}
	return alfa;
}