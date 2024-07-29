#include "Basis.h"
#include "AlgebraHelper.h"

double Basis::BasisFunction::Spline::GetValue(double x) {
	return a + b * x + c * x * x + d * x * x * x;
}

double Basis::BasisFunction::Spline::GetDerivative(double x) {
	return b + 2.0 * c * x + 3.0 * d * x * x;
}

Basis::BasisFunction::Spline::Spline(double A, double B, double C, double D) {
	this->a = A;
	this->b = B;
	this->c = C;
	this->d = D;
}

void Basis::BasisFunction::AddSpline(double A, double B, double C, double D) {
	_splines.push_back(Spline(A, B, C, D));
}

double Basis::BasisFunction::GetDerivative(int spline, double relativeDistance) {
	return _splines[spline].GetDerivative(relativeDistance);
}

double Basis::BasisFunction::GetValue(int spline, double relativeDistance) {
	return _splines[spline].GetValue(relativeDistance);
}

Basis::Basis(int Points) {
	_points = Points;

	auto h = std::make_unique<double[]>(_points - 1);
	for (int i = 0; i < _points - 1; ++i) {
		h[i] = 1.0;
	}

	auto ah = std::make_unique<AlgebraHelper>();
	auto M = ah->GenerateTriDiagonal(_points, h);

	auto R = std::make_unique<std::unique_ptr<double[]>[]>(_points);
	for (int k = 0; k < _points; ++k) {
		R[k] = std::make_unique<double[]>(_points);
	}
	ah->gaussJordan(M, R, _points);

	for (int i = 0; i < _points; ++i) {
		auto e = std::make_unique<double[]>(_points);
		for (int j = 0; j < _points; ++j) {
			e[j] = 0.0;
		}
		e[i] = 1.0;

		int len = _points - 1;
		auto a = std::make_unique<double[]>(len);
		auto b = std::make_unique<double[]>(len);
		auto c = std::make_unique<double[]>(len);
		auto d = std::make_unique<double[]>(len);

		ah->MakeSplines(R, e, h, a, b, c, d, _points);

		BasisFunction basicFunction;
		for (int j = 0; j < len; ++j) {
			basicFunction.AddSpline(a[j], b[j], c[j], d[j]);
		}
		_basisFunctions.push_back(basicFunction);
	}
}

int Basis::nPoints() {
	return _points;
}

std::unique_ptr<double[]> Basis::GetAllValues(int k, double relative) {
	auto values = std::make_unique<double[]>(_points);
	for (int i = 0; i < _points; ++i) {
		values[i] = _basisFunctions[i].GetValue(k, relative);
	}
	return values;
}

std::unique_ptr<double[]> Basis::GetAllDerivatives(int k, double relative) {
	auto derivative = std::make_unique<double[]>(_points);
	for (int i = 0; i < _points; ++i) {
		derivative[i] = _basisFunctions[i].GetDerivative(k, relative);
	}
	return derivative;
}
