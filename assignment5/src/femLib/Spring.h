#pragma once

#include "Element.h"

/**
	This class implements the interface for an elementary energy unit. As a function of deformed, undeformed,
	and other parameters, such as boundary conditions, each class that extends this one will define a potential energy.
	The deformed energy depends on a number of nodes.
*/
class Spring : public Element {

public:
	Spring(const std::array<int, 2> &nodeIndices, const VectorXd &X)
		: nodeIndices(nodeIndices) {
	}
	virtual ~Spring() {}

	// Returns the number of nodes this unit depends on
	virtual int getNumNodes() const {
		return 2;
	}
	// Returns the global index of node `i`
	virtual int getNodeIndex(int i) const {
		return nodeIndices[i];
	}

	// Returns the element's mass
	virtual double getMass() const {
		return 0;
	}

	// Returns the energy value given deformed `x` and undeformed `X` state
	virtual double getEnergy(const VectorXd& x, const VectorXd& X) {

		// Ex 1.2
		// Task: Given `x` and `X`, return the spring energy.

		// Some notes:
		// `x` and `X` contain the current and rest positions of all
		// nodes. You can extract the position of e.g. node 0 like this:
		// Vector2d x1 = getVertex(0, x);
		// or to get the rest position of node 0:
		// Vector X1 = getVertex(0, X);
		// The spring stiffness is stored in the variable `k`.

		double l = (getNodePos(1, x) - getNodePos(0, x)).norm();

		double L = (getNodePos(1, X) - getNodePos(0, X)).norm();

		double energy = (k * pow((l / L) - 1, 2) * L) / 2;

		return energy;

	}

	// Adds the gradient to `grad` given deformed `x` and undeformed `X` state
	virtual void addEnergyGradientTo(const VectorXd& x, const VectorXd& X, VectorXd& grad) {

		// Ex 1.2
		// Task: Given `x` and `X`, add the gradient of the spring energy to `grad`.

		// Again, you can extract the position of e.g. node 0 like this:
		// Vector2d x1 = getVertex(0, x);
		// and the spring stiffness is stored in `k`.

		// Remember that `grad` is a vector of size 2*N, where N is the total
		// number of nodes in the system. Make sure you are writing to the
		// correct location in `grad`. To get the global index of node 0 of
		// this spring, use this function:
		// int globalIndex0 = getNodeIndex(0);
		// or for node 1
		// int globalIndex1 = getNodeIndex(1);

		// Original and new lengths

		double L = (getNodePos(0, X) - getNodePos(1, X)).norm();

		double l = (getNodePos(0, x) - getNodePos(1, x)).norm();

		// Force

		Vector2d f = k * (1 - l / L) * (getNodePos(1, x) - getNodePos(0, x)) / l;

		// Gradients

		grad[2 * getNodeIndex(0)] += f[0];
		grad[2 * getNodeIndex(0) + 1] += f[1];

		grad[2 * getNodeIndex(1)] += -f[0];
		grad[2 * getNodeIndex(1) + 1] += -f[1];

	}

	// Adds the hessian entries to `hesEntries` given deformed `x` and undeformed `X` state
	virtual void addEnergyHessianTo(const VectorXd& x, const VectorXd& X, std::vector<Tripletd>& hesEntries) {

		// Ex 1.4
		// Task: Given `x` and `X`, add the hessian of the spring energy to `hesEntries`.

		// Original and new lengths

		double L = (getNodePos(0, X) - getNodePos(1, X)).norm();

		double l = (getNodePos(0, x) - getNodePos(1, x)).norm();

		Vector2d l_vector = getNodePos(1, x) - getNodePos(0, x);

		// Second derivative

		double dE_dx1dx1 = (k * (l / L - 1) / l) * (1 - pow(l_vector(0) / l, 2)) + k * pow(l_vector(0) / l, 2) / L;

		double dE_dx1dy1 = l_vector(0) / l * l_vector(1) / l * ((k / l) - (k * (l / L - 1) / l));

		double dE_dy1dy1 = (k * (l / L - 1) / l) + ((k / l) - (k * (l / L - 1) / l)) * pow(l_vector(1) / l, 2);

		// Gradients

		double dE_dx1dx2 = -dE_dx1dx1;
		double dE_dx1dy2 = -dE_dx1dy1;

		double dE_dy1dx1 = dE_dx1dy1;
		double dE_dy1dx2 = -dE_dx1dy1;
		double dE_dy1dy2 = -dE_dy1dy1;

		double dE_dx2dx1 = -dE_dx1dx1;
		double dE_dx2dy1 = -dE_dx1dy1;
		double dE_dx2dx2 = dE_dx1dx1;
		double dE_dx2dy2 = dE_dx1dy1;

		double dE_dy2dx1 = -dE_dx1dy1;
		double dE_dy2dy1 = -dE_dy1dy1;
		double dE_dy2dx2 = dE_dx1dy1;
		double dE_dy2dy2 = dE_dy1dy1;

		// Hessian entries

		double nodeIndex_x1 = 2 * getNodeIndex(0);
		double nodeIndex_y1 = nodeIndex_x1 + 1;
		double nodeIndex_x2 = 2 * getNodeIndex(1);
		double nodeIndex_y2 = nodeIndex_x2 + 1;

		hesEntries.push_back(Tripletd(nodeIndex_x1,		nodeIndex_x1,		dE_dx1dx1));
		hesEntries.push_back(Tripletd(nodeIndex_x1,		nodeIndex_y1,		dE_dx1dy1));
		hesEntries.push_back(Tripletd(nodeIndex_x1,		nodeIndex_x2,		dE_dx1dx2));
		hesEntries.push_back(Tripletd(nodeIndex_x1,		nodeIndex_y2,		dE_dx1dy2));

		hesEntries.push_back(Tripletd(nodeIndex_y1,		nodeIndex_x1,		dE_dy1dx1));
		hesEntries.push_back(Tripletd(nodeIndex_y1,		nodeIndex_y1,		dE_dy1dy1));
		hesEntries.push_back(Tripletd(nodeIndex_y1,		nodeIndex_x2,		dE_dy1dx2));
		hesEntries.push_back(Tripletd(nodeIndex_y1,		nodeIndex_y2,		dE_dy1dy2));

		hesEntries.push_back(Tripletd(nodeIndex_x2,		nodeIndex_x1,		dE_dx2dx1));
		hesEntries.push_back(Tripletd(nodeIndex_x2,		nodeIndex_y1,		dE_dx2dy1));
		hesEntries.push_back(Tripletd(nodeIndex_x2,		nodeIndex_x2,		dE_dx2dx2));
		hesEntries.push_back(Tripletd(nodeIndex_x2,		nodeIndex_y2,		dE_dx2dy2));

		hesEntries.push_back(Tripletd(nodeIndex_y2,		nodeIndex_x1,		dE_dy2dx1));
		hesEntries.push_back(Tripletd(nodeIndex_y2,		nodeIndex_y1,		dE_dy2dy1));
		hesEntries.push_back(Tripletd(nodeIndex_y2,		nodeIndex_x2,		dE_dy2dx2));
		hesEntries.push_back(Tripletd(nodeIndex_y2,		nodeIndex_y2,		dE_dy2dy2));

	}

protected:
	// the collection of nodes that define the triangle element
	std::array<int, 2> nodeIndices;
	// spring stiffness
	double k = 20.0;
};
