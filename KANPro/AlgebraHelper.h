#pragma once
#include <memory>

class AlgebraHelper
{
private:
	bool choldc1(int n, std::unique_ptr<std::unique_ptr<double[]>[]>& a, std::unique_ptr<double[]>& p);
	bool choldc(int n, std::unique_ptr<std::unique_ptr<double[]>[]>& A, std::unique_ptr<std::unique_ptr<double[]>[]>& a);
	void ShowMatrix(std::unique_ptr<std::unique_ptr<double[]>[]>& matrix, int rows, int cols);
	void MultiplySquare(std::unique_ptr<std::unique_ptr<double[]>[]>& X, std::unique_ptr<std::unique_ptr<double[]>[]>& Y,
		std::unique_ptr<std::unique_ptr<double[]>[]>& Z, int size);
public:
	bool CholeskySolution(std::unique_ptr<std::unique_ptr<double[]>[]>& M, std::unique_ptr<std::unique_ptr<double[]>[]>& Inv, int size);
	void InverseSelfTest();
	std::unique_ptr<std::unique_ptr<double[]>[]> GenerateTriDiagonal(int N, std::unique_ptr<double[]>& h);
	void SplinesSelfTest();
	bool gaussJordan(std::unique_ptr<std::unique_ptr<double[]>[]>& matrix, std::unique_ptr<std::unique_ptr<double[]>[]>& inverse, int size);
	void MakeSplines(std::unique_ptr<std::unique_ptr<double[]>[]>& A,
		std::unique_ptr<double[]>& y, std::unique_ptr<double[]>& h,
		std::unique_ptr<double[]>& a, std::unique_ptr<double[]>& b, std::unique_ptr<double[]>& c, std::unique_ptr<double[]>& d, int lenY);
};

