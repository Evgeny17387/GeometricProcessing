#include <iostream>

#include "optLib/GradientDescentFunctionMinimizer.h"
#include "optLib/NewtonFunctionMinimizer.h"
#include "optLib/RosenbrockFunction.h"

int main(int argc, char *argv[]) {

	RosenbrockFunction function;

	VectorXd x(2); x << 1.1, 1.1;

#if 1
	// Ex 1.1: Gradient descent
	{
		GradientDescentFunctionMinimizer minimizer(10000);

		bool ok = minimizer.minimize(&function, x);

		std::cout << "---\nGradient descent:" << std::endl;
		std::cout << "converged:  " << ((ok) ? "yes" : "no") << std::endl;
		std::cout << "iterations: " << minimizer.getLastIterations() << std::endl;
		std::cout << "min. x:     " << x.transpose() << std::endl;
		std::cout << "min. value: "<< function.computeValue(x) << std::endl;
	}
#endif

	// Ex 1.3: Newton's method
#if 0 // set to 1 for Ex 1.3
	{
		NewtonFunctionMinimizer minimizer;

		bool ok = minimizer.minimize(&function, x);

		std::cout << "---\nNewton's method:" << std::endl;
		std::cout << "converged:  " << ((ok) ? "yes" : "no") << std::endl;
		std::cout << "iterations: " << minimizer.getLastIterations() << std::endl;
		std::cout << "min. x:     " << x.transpose() << std::endl;
		std::cout << "min. value: "<< function.computeValue(x) << std::endl;
	}
#endif

}
