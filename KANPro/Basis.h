#pragma once
#include <vector>
#include <memory>

struct Basis {
	struct BasisFunction {
		struct Spline {
			double a, b, c, d;
			double GetValue(double x);
			double GetDerivative(double x);
			Spline(double A, double B, double C, double D);
		};
		std::vector<Spline> _splines;
		void AddSpline(double A, double B, double C, double D);
		double GetDerivative(int spline, double relativeDistance);
		double GetValue(int spline, double relativeDistance);
	};

	int _points = 0;
	std::vector<BasisFunction> _basisFunctions;
	Basis(int Points);
	int nPoints();
	std::unique_ptr<double[]> GetAllValues(int k, double relative);
	std::unique_ptr<double[]> GetAllDerivatives(int k, double relative);
};