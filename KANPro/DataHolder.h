#pragma once
#include <iostream>

class DataHolder {
public:
	bool ReadDataWDBC();
	bool ReadDataSpam();
	bool ReadDataMFormula();
	int nSeparatedBlocks(std::string data, char delimiter);
	std::unique_ptr<std::unique_ptr<double[]>[]> inputs;
	std::unique_ptr<double[]> target;
	int nRecords = -1;
	int nFeatures = -1;
};

