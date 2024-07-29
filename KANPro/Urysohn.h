#pragma once
#include <memory>
#include <vector>
#include "Univariate.h"

class Urysohn {
private:
	int _nPoints;
	std::vector<std::unique_ptr<Univariate>> _univariateList;
public:
	Urysohn(std::unique_ptr<double[]>& xmin, std::unique_ptr<double[]>& xmax, double targetMin, 
		double targetMax, Basis& basis, int nFeatures); 
	double GetDerivative(int layer, double x);
	void UpdateUsingMemory(double delta, double mu);
	void UpdateUsingInput(double delta, std::unique_ptr<double[]>& inputs, double mu);
	double GetValueUsingInput(std::unique_ptr<double[]>& inputs);
	void AssignUPoints(std::unique_ptr<double[]>& c, int n, int len);
};

