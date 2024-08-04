#include"utilities.h"

//Basic legendre polynomials

long double Legendre(unsigned short int poly_degree, long double x) {
	if (poly_degree == 0) {
		return 1.0;
	}
	else if (poly_degree == 1) {
		return x;
	}
	else {
		long double Pn_1 = 1.0; // P_0(x)
		long double Pn = x;     // P_1(x)
		long double temp;

		for (unsigned int i = 1; i < poly_degree; ++i) {
			temp = ((2.0 * i + 1.0) * x * Pn - i * Pn_1) / (i + 1.0);
			Pn_1 = Pn;
			Pn = temp;
		}
		return Pn;
	}
}

//Derivative of the legendre polynomials
// 
 //from the formula: dP=n*((Pn_1-x*Pn)/(1-x^2))
			 /*long double Pn, Pn_1;
			 Legendre_d_helper(poly_degree, x, Pn, Pn_1);
			 return poly_degree * (Pn_1 - x * Pn) / (1.0L - x * x);*/

 long double Legendre_d(unsigned short int poly_degree, long double x) {
	 try {
		 if (poly_degree < 1) {
			 throw std::runtime_error("Invalid poly degree, there is no 0th degree derivative poly here!");
		 }
		 else {
			 switch (poly_degree) {
			 case 1:
				 return 1.0L;
				 break;
			 case 2:
				 return 3.0L * x;
				 break;
			 default:
				 long double dP_1 = 3.0*x;
				 long double dP = 0.0L;
				 long double temp = 0.0L;
				 for (unsigned short int i = 3; i < poly_degree+1; i++) {
					 temp = i * Legendre(i-1, x) + x * dP_1;
					 dP = temp;
					 dP_1 = dP;
				 }
				 return dP;
				 break;
			 }
		 }
	 }
	 catch (const std::runtime_error& error) {
		 std::cerr << "Runtime error in Legendre_d:" << error.what() << std::endl;
		 return 0;
	 }
}


//integrated legendre polynomials

 long double integrated_Leg(unsigned short int poly_degree, long double x) {
	switch (poly_degree)
	{
	case 1:
		return(0.5L * (1.0L - x));
		break;
	case 2:
		return(0.5L * (1.0L + x));
		break;
	case 3: 
		return (3.0L / (2.0L * std::sqrt(6.0L))) * (x * x - 1.0L);
		break;
	default:
		return (1.0L / std::sqrt(4.0L * (poly_degree - 1) - 2.0L)) * (Legendre(poly_degree - 1, x) - Legendre(poly_degree - 3, x));
		break;
	}
}

//derivative of the integrated legendre polynomials

 long double integrated_Leg_d(unsigned short int poly_degree, long double x) {
	//Legendre_d handles runtime error here
	 switch (poly_degree)
	{
	case 1:
		return(-0.5L);
		break;
	case 2:
		return(0.5L);
		break;
	case 3:
		return(3.0L/std::sqrt(6.0L)*x);
		break;
	default:
		long double factor = 1.0L / std::sqrt(4.0L * (poly_degree - 1) - 2.0L);
		return factor * (Legendre_d(poly_degree - 1, x) - Legendre_d(poly_degree - 3, x));		
		break;
	}
}

//partitioning an interval into "pieces" pieces, result is in a vector

std::vector<long double> partitioner(long double interval_start, long double interval_end, unsigned int piece) { 
	if (interval_start > interval_end) {
		std::swap(interval_start, interval_end);
	}
	long double h= 0.0L;
	h = std::abs(interval_end - interval_start) / piece;
	std::vector<long double> result(piece+1);
	result[0]=interval_start;
	for (unsigned int i = 1; i <= piece; i++) {
		result[i] = interval_start + i * h;
	}
	//result[piece] = interval_end;
	return result;
}

//Exponential heating test function

Eigen::VectorXd Test_function(long double time, unsigned short int poly_degree, unsigned int n) {

	if (time < 0) {
		throw std::runtime_error("There is no negative time, you sunova...!"); //doin' it for fun and practice
	}
	else if (poly_degree < 1 || n < 1) {
		throw std::runtime_error("Polydegree or element number is not enough, try again!");
	}
	else {
		Eigen::VectorXd F = Eigen::VectorXd::Zero(2 * (poly_degree * n - n + 1));
		F(0) = (-10000 * (1 / 0.075) * 6 / (6 - (1 / 0.075)) * (exp((-1 / 0.075) * time / 0.008) - exp(-6 * time / 0.008)));
		return F;
	}	
}

//Laser heating test function

Eigen::VectorXd Laser_heating(long double time, unsigned short int poly_degree, unsigned int n) {
	Eigen::VectorXd F = Eigen::VectorXd::Zero(2 * (poly_degree * n - n + 1));
	if ((0.001 < time && time < 2.5) || (5 < time && time < 7.5)) {
		F(0) = -140000;
	}
	else {
		F(0) = 0;
	}
	return F;
}

// f vector for Test_function

Eigen::VectorXd f_test(unsigned short int poly_degree, unsigned int n, double theta, long double time, long double dt) {
	
	//The function input could be any  test function with the desired parameters

	//Declaration of eigen::vectors
	Eigen::VectorXd F = Test_function(time, poly_degree, n);
	Eigen::VectorXd dF = Test_function(time+dt, poly_degree, n);
	Eigen::VectorXd f = Eigen::RowVectorXd::Zero(2 * (poly_degree * n - n + 1));

	f(0) = (F(0) * (1 - theta) + dF(0) * theta) * dt;
	return f;
} 

//Calculate nodes and weights for the gaussian quadrature, works perfectly
Eigen::MatrixXd gauss_legendre(unsigned short int poly_degree) {
	const long double PI = 3.14159265358979323846264338L;
		long double p1=0.0, p2=0.0, p3 = 0.0, pp = 0.0, z1 = 0.0;
		Eigen::MatrixXd nodesandWeights=Eigen::MatrixXd::Zero(poly_degree,2);
	for (int i = 0; i < poly_degree; ++i) {
		z1 = cos(PI * (i + 0.75) / (poly_degree + 0.5));
		do {
			p1 = 1.0;
			p2 = 0.0;
			for (int j = 0; j < poly_degree; ++j) {
				p3 = p2;
				p2 = p1;
				p1 = ((2.0 * j + 1.0) * z1 * p2 - j * p3) / (j + 1);
			}
			pp = poly_degree * (z1 * p1 - p2) / (z1 * z1 - 1.0);
			z1 -= p1 / pp;
		} while (fabs(p1 / pp) > 1e-15);
		nodesandWeights(i,0) = z1;
		nodesandWeights(i,1) = 2.0 / ((1.0 - z1 * z1) * pp * pp);
	}
	return nodesandWeights;
}

long double GaussianQuadrature(std::function<long double(unsigned short int, long double)>func1, std::function<long double(unsigned short int, long double)>func2,
	unsigned short int poly_degree1, unsigned short int poly_degree2) {
	unsigned short int total_degree = (poly_degree1 + poly_degree2);
	long double sum = 0.0L;
	Eigen::MatrixXd NodesandWeights = gauss_legendre(total_degree+1);
	for (unsigned short int i = 0; i < total_degree+1; i++) {
		sum += NodesandWeights(i, 1) * (func1(poly_degree1, NodesandWeights(i, 0)) * func2(poly_degree2, NodesandWeights(i, 0)));
	}
	return sum;
}