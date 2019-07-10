#pragma once

#include "ObjectiveFunction.h"

#include <string>
#include <iostream>
#include <cmath>
#include <cfloat>

class GradientDescentFunctionMinimizer{
public:
	GradientDescentFunctionMinimizer(int maxIterations=100, double solveResidual=1e-5, int maxLineSearchIterations=15)
		: maxIterations(maxIterations), solveResidual(solveResidual), maxLineSearchIterations(maxLineSearchIterations){
	}

	virtual ~GradientDescentFunctionMinimizer(){}

	int getLastIterations() { return lastIterations; }

	virtual bool minimize(ObjectiveFunction *function, VectorXd &x){

		//number of parameters...
		int N = (int) x.size();
		resize(xi, N);
		resize(dx, N);
		resize(gradient, N);

		xi = x;

		bool optimizationConverged = false;

		int i=0;
		for(; i < maxIterations; i++) {

			std::cout << "Iteration: " << i << std::endl;

			computeSearchDirection(function, xi, dx);

			if (dx.norm() < solveResidual){
				optimizationConverged = true;
				break;
			}

			doLineSearch(function, dx, xi);
		}

		lastIterations = i;

		//p now holds the parameter values at the start of the iteration...
		x = xi;

		//and done!
		return optimizationConverged;
	}

protected:
	// Since the gradient of a function gives the direction of steepest descent, all one needs to do is go in that direction...
	virtual void computeSearchDirection(ObjectiveFunction *function, const VectorXd &x, VectorXd& dx) {

		// Ex. 1.1

		// Init step

		resize(dx, x.size());

		// Assumes step = 1

		function->addGradientTo(dx, x);
		dx = -dx;

	}

	virtual void doLineSearch(ObjectiveFunction *function, const VectorXd& dx, VectorXd& xi)
	{

		// Ex. 1.1

		// Init

		double value = function->computeValue(xi);

		double alphaCurrent = alpha;

		double valueNew = function->computeValue(xi + alphaCurrent * dx);
		
		// Loop until finds right step or exceeds mac number of iterations for line search

		int i = 0;
		while ( !(valueNew < value) && (i < maxLineSearchIterations) )
		{

			alphaCurrent = alphaCurrent * beta;

			valueNew = function->computeValue(xi + alphaCurrent * dx);

			i++;

		}

		xi = xi + alphaCurrent * dx;

	}

protected:

	double solveResidual = 1e-5;
	int maxIterations = 100;
	int maxLineSearchIterations = 15;

	VectorXd xi, dx, gradient;

	// some stats about the last time `minimize` was called
	int lastIterations = -1;

	double alpha = 1.0;
	double beta = 0.5;

};
