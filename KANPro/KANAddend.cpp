#include "KANAddend.h"

KANAddend::KANAddend(std::unique_ptr<double[]>& xmin, std::unique_ptr<double[]>& xmax, double targetMin, double targetMax,
	Basis& innerBasis, Basis& outerBasis, double mu_inner, double mu_outer, int nFeatures) {
    _muInner = mu_inner;
    _muOuter = mu_outer;
    _univariate = std::make_unique<Univariate>(targetMin, targetMax, targetMin, targetMax, outerBasis);
    _u = std::make_unique<Urysohn>(xmin, xmax, targetMin, targetMax, innerBasis, nFeatures);
}

void KANAddend::UpdateOuterPoints(std::unique_ptr<double[]>& c, int len) {
    _univariate->AssignCoefficients(c, len);
}

void KANAddend::UpdateUsingMemory(double diff) {
    double derivative = _univariate->GetDerivative(_lastInnerValue);
    _u->UpdateUsingMemory(diff * derivative, _muInner);
    _univariate->UpdateUsingMemory(diff, _muOuter);
}

void KANAddend::UpdateUsingInput(std::unique_ptr<double[]>& inputs, double diff) {
    double value = _u->GetValueUsingInput(inputs);
    double derivative = _univariate->GetDerivative(value);
    _u->UpdateUsingInput(diff * derivative, inputs, _muInner);
    _univariate->UpdateUsingInput(value, diff, _muOuter);
}

double KANAddend::ComputeUsingInput(std::unique_ptr<double[]>& inputs) {
    _lastInnerValue = _u->GetValueUsingInput(inputs);
    return _univariate->GetFunctionUsingInput(_lastInnerValue);
}
