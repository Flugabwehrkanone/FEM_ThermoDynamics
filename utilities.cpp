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

//Performace booster for Legendre_d function

//void Legendre_d_helper(unsigned short int poly_degree, long double x, long double& Pn, long double& Pn_1) {
//	if (poly_degree == 0) {
//		Pn = 1.0L;
//		Pn_1 = 0.0L; // P_{-1} is not defined, but this value will not be used
//	}
//	else if (poly_degree == 1) {
//		Pn = x;
//		Pn_1 = 1.0L; // P_0(x) = 1
//	}
//	else {
//		Pn_1 = 1.0L; // P_0(x)
//		Pn = x;     // P_1(x)
//		long double Pn_minus_2 = Pn_1; // To hold P_{n-2} during iteration
//
//		for (unsigned short int i = 1; i < poly_degree; ++i) {
//			long double Pn_temp = Pn;
//			Pn = ((2.0L * i + 1.0L) * x * Pn - i * Pn_minus_2) / (i + 1.0L);
//			Pn_minus_2 = Pn_temp;
//			Pn_1 = Pn_temp;
//		}
//	}
//}

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

//void swapColumns(std::vector<std::vector<long double>>& matrix, unsigned short int col1, unsigned short int col2) {
//	for (unsigned int i = 0; i < matrix.size(); ++i) {
//		std::swap(matrix[i][col1], matrix[i][col2]);
//	}
//}
//
//void Matrix_multiply(std::vector<std::vector<long double>>& matrix, float multiplier) {
//	for (auto& row : matrix) {
//		for (auto& elem : row) {
//			elem *= multiplier;
//		}
//	}
//}
// 
 // Simpson's Rule algorithm
 
//Simpson's rule, works perfectly
long double Simpsons_rule(const std::vector<long double>& functionValues, long double a, long double b) {
	 // Check if n is even
	unsigned int n = functionValues.size()-1;
	 if (n % 2 != 0) {
		 throw std::invalid_argument("Number of intervals n must be even.");
	 }

	 // Compute the width of each interval
	 double h = (b - a) / n;
	 double integral = functionValues[0] + functionValues[n];

	 // Apply Simpson's rule
#pragma omp parallel for
	 for (int i = 1; i < n; i += 2) {
		 integral += 4 * functionValues[i];
	 }
#pragma omp parallel for
	 for (int i = 2; i < n - 1; i += 2) {
		 integral += 2 * functionValues[i];
	 }
	 integral *= h / 3;
	 return integral;
 }

//makes an interval with the functions values for the simpsons rule algorithm
//the algorithm assumes that the difference between the intervals start and end is 2
//this is so, that the function can have less arguments
//this function works perfectly

std::vector<long double>Function_integr_values(std::function<long double(unsigned short int, long double)> func1, std::function<long double(unsigned short int, long double)> func2, 
unsigned int poly_degree1, unsigned int poly_degree2) {
	unsigned int number_of_points = ceil((poly_degree1 * poly_degree2) / 2.0);
	
	if (number_of_points < 3) {
		number_of_points = 3;
	}
	if (number_of_points % 2 == 0) {
		number_of_points += 1;
	}
	
	std::vector<long double>interval(number_of_points,0);
	long double h = 2.0L / (number_of_points-1);
	long double x = 0.0;
	
#pragma omp parallel for
	for (int i = 0; i < number_of_points; i++) {
		x = -1.0L + i * h;
		//std::cout << x << std::endl;
		interval[i] = func1(poly_degree1, x) * func2(poly_degree2, x);
	}
	return interval;
}