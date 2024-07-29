#include "Urysohn.h"

Urysohn::Urysohn(std::unique_ptr<double[]>& xmin, std::unique_ptr<double[]>& xmax, double targetMin, 
	double targetMax, Basis& basis, int nFeatures) {
	_nPoints = nFeatures;
	double ymin = targetMin / _nPoints;
	double ymax = targetMax / _nPoints;
	for (int i = 0; i < _nPoints; ++i)
	{
		_univariateList.push_back(std::make_unique<Univariate>(xmin[i], xmax[i], ymin, ymax, basis));
	}
}

double Urysohn::GetDerivative(int layer, double x) {
	return _univariateList[layer].get()->GetDerivative(x);
}

void Urysohn::UpdateUsingMemory(double delta, double mu) {
	for (int i = 0; i < _nPoints; ++i) {
		_univariateList[i].get()->UpdateUsingMemory(delta, mu);
	}
}

void Urysohn::UpdateUsingInput(double delta, std::unique_ptr<double[]>& inputs, double mu) {
	for (int i = 0; i < _nPoints; ++i) {
		_univariateList[i].get()->UpdateUsingInput(inputs[i], delta, mu);
	}
}
	
double Urysohn::GetValueUsingInput(std::unique_ptr<double[]>& inputs) {
	double f = 0.0;
	for (int i = 0; i < _nPoints; ++i) {
		f += _univariateList[i].get()->GetFunctionUsingInput(inputs[i]);
	}
	return f;
}

void Urysohn::AssignUPoints(std::unique_ptr<double[]>& c, int n, int len) {
	_univariateList[n]->AssignCoefficients(c, len);
}