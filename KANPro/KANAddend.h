#pragma once
#include "Urysohn.h"

class KANAddend
{
public:
	KANAddend(std::unique_ptr<double[]>& xmin, std::unique_ptr<double[]>& xmax, double targetMin, double targetMax,
		Basis& innerBasis, Basis& outerBasis, double mu_inner, double mu_outer, int nFeatures);
	void UpdateUsingMemory(double diff);
	void UpdateUsingInput(std::unique_ptr<double[]>& inputs, double diff);
	double ComputeUsingInput(std::unique_ptr<double[]>& inputs);
	void UpdateOuterPoints(std::unique_ptr<double[]>& c, int len);
	std::unique_ptr<Urysohn> _u;
private:
	std::unique_ptr<Univariate> _univariate;
	double _muInner, _muOuter, _lastInnerValue;
};

