	#pragma once

#include "ObjectiveFunction.h"
#include "GradientDescentFunctionMinimizer.h"

/**
	use Newton's method to optimize a function. p will store the final value that minimizes the function, and its initial value
	is used to start the optimization method.

	Task: find p that minimize f(p). This means that df/dp(p) = 0.
	df/dp(p+dp) ~ df/dp(p) + d/dp(df/dp) * dp = 0 ==> -df/dp(p) = d/dp(df/dp) * dp
	Iterating the above, will hopefully get p that minimizes f.
*/
class NewtonFunctionMinimizer : public GradientDescentFunctionMinimizer {
public:
	NewtonFunctionMinimizer(int maxIterations = 100, double solveResidual = 0.0001, int maxLineSearchIterations = 15)
		: GradientDescentFunctionMinimizer(maxIterations, solveResidual, maxLineSearchIterations) {	}

	virtual ~NewtonFunctionMinimizer() {}

protected:
	// The search direction is given by -Hinv * g
	virtual void computeSearchDirection(ObjectiveFunction *function, const VectorXd &x, VectorXd& dx) {

		// Ex 1.3

		// Compute Hessian

		std::vector<Tripletd> hessianEntries;
		function->addHessianEntriesTo(hessianEntries, x);
		Eigen::SparseMatrix<double> hessian(x.rows(), x.rows());
		hessian.setFromTriplets(hessianEntries.begin(), hessianEntries.end());

		// Compute Gradient

		VectorXd gradient;
		function->addGradientTo(gradient, x);
		gradient = -gradient;

		// Solve equation

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		solver.compute(hessian);

		dx = solver.solve(gradient);

	}

public:
	SparseMatrixd H;
	std::vector<Tripletd> hessianEntries;
};
