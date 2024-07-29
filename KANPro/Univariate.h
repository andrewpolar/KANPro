#pragma once
#include <memory>
#include "Basis.h"

class Univariate
{
private:
	int _points;
	double _xmin, _xmax, _ymin, _ymax, _deltax;
	std::unique_ptr<double[]> _vlast; 
	std::unique_ptr<double[]> _coefficients;
	Basis& _basis;
	void SetLimits();
	void FitDefinition(double x);
public:
	Univariate(double xmin, double xmax, double ymin, double ymax, Basis& basis);
	void GetSplineIndexAndOffset(double x, int& k, double& relative);
	double GetDerivative(double x);
	double GetFunctionUsingInput(double x);
	void UpdateUsingMemory(double delta, double mu);
	void UpdateUsingInput(double x, double delta, double mu);
	void AssignCoefficients(std::unique_ptr<double[]>& c, int len);
};

