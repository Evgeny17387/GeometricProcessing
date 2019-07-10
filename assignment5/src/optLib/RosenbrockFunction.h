#pragma once

#include "ObjectiveFunction.h"

class RosenbrockFunction : public ObjectiveFunction {
public:

    RosenbrockFunction() {
		a = 1;
		b = 100;
    }

    virtual double computeValue(const VectorXd& x) {

		// Ex 1.1

		return pow((a - x(0)), 2) + b * pow(x(1) - pow(x(0), 2), 2);

	}

    virtual void addGradientTo(VectorXd& grad, const VectorXd& x) {

		// Ex 1.1

		resize(grad, x.size());

		grad(0) = -2 * (a - x(0)) - 4 * b * x(0) * (x(1) - pow(x(0), 2));
		grad(1) = 2 * b * (x(1) - pow(x(0), 2));

	}

	virtual void addHessianEntriesTo(std::vector<Tripletd>& hessianEntries, const VectorXd& x) {

		// Ex 1.2

		hessianEntries.push_back(Tripletd(0, 0, 2 - 4 * b * x(1) + 12 * b * pow(x(0), 2)));
		hessianEntries.push_back(Tripletd(0, 1, -4 * b * x(0)));
		hessianEntries.push_back(Tripletd(1, 0, -4 * b * x(0)));
		hessianEntries.push_back(Tripletd(1, 1, 2 * b));

	}

    double a, b;

};
