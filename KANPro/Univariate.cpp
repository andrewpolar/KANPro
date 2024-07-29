#include "Univariate.h"

void Univariate::SetLimits() {
	double range = _xmax - _xmin;
	_xmin -= 0.01 * range;
	_xmax += 0.01 * range;
	_deltax = (_xmax - _xmin) / (_points - 1);
}
	
void Univariate::FitDefinition(double x) {
	if (x < _xmin) {
		_xmin = x;
		SetLimits();
	}
	if (x > _xmax) {
		_xmax = x;
		SetLimits();
	}
}

Univariate::Univariate(double xmin, double xmax, double ymin, double ymax, Basis& basis) : _basis(basis) {
	_points = _basis.nPoints();
	_xmin = xmin;
	_xmax = xmax;
	SetLimits();
	_ymin = ymin;
	_ymax = ymax;
	_coefficients = std::make_unique<double[]>(_points);
	for (int i = 0; i < _points; i++) {
		_coefficients[i] = (rand() % 990 + 10) / 1000.0 * (_ymax - _ymin) + _ymax;
		_coefficients[i] /= _points;
	}
}

void Univariate::GetSplineIndexAndOffset(double x, int& k, double& relative) {
	double offset = (x - _xmin) / _deltax;
	k = (int)(offset);
	relative = offset - k;
}
	
double Univariate::GetDerivative(double x) {
	FitDefinition(x);
	int k; 
	double relative;
	GetSplineIndexAndOffset(x, k, relative);
	auto derivative = _basis.GetAllDerivatives(k, relative);
	double v = 0.0;
	for (int i = 0; i < _points; i++) {
		v += derivative[i] * _coefficients[i];
	}
	return v / _deltax;
}
	
double Univariate::GetFunctionUsingInput(double x) {
	FitDefinition(x);
	int k;
	double relative;
	GetSplineIndexAndOffset(x, k, relative);
	_vlast = _basis.GetAllValues(k, relative);
	double f = 0.0;
	for (int i = 0; i < _points; i++) {
		f += _vlast[i] * _coefficients[i];
	}
	return f;
}

void Univariate::UpdateUsingMemory(double delta, double mu) {
	for (int i = 0; i < _points; ++i) {
		_coefficients[i] += delta * mu * _vlast[i];
	}
}

void Univariate::UpdateUsingInput(double x, double delta, double mu) {
	FitDefinition(x);
	int k;
	double relative;
	GetSplineIndexAndOffset(x, k, relative);
	_vlast = _basis.GetAllValues(k, relative);
	for (int i = 0; i < _points; ++i) {
		_coefficients[i] += delta * mu * _vlast[i];
	}
}

void Univariate::AssignCoefficients(std::unique_ptr<double[]>& c, int len) {
	for (int i = 0; i < len; ++i) {
		_coefficients[i] = c[i];
	}
}
