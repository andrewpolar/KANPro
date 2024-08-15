#include "DataHolder.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

int DataHolder::nSeparatedBlocks(std::string data, char delimiter) {
	std::string::difference_type n = std::count(data.begin(), data.end(), delimiter);
	return (int)n;
}

bool DataHolder::ReadDataWDBC() {
	std::ifstream myfile;
	myfile.open("..//DATA//WDBC//wdbc.data");
	std::string delimiter = ",";
	nRecords = 570;
	nFeatures = 30;

	target = std::make_unique<double[]>(nRecords);
	inputs = std::make_unique<std::unique_ptr<double[]>[]>(nRecords);
	for (int k = 0; k < nRecords; ++k) {
		inputs[k] = std::make_unique<double[]>(nFeatures);
	}

	std::string myline;
	auto data = std::make_unique<std::string[]>(32);
	int nLine = 0;
	if (myfile.is_open()) {
		while (myfile.good()) {
			myfile >> myline;
			myline += ",";
			size_t pos = 0;
			std::string token;
			int i = 0;
			while ((pos = myline.find(delimiter)) != std::string::npos) {
				token = myline.substr(0, pos);
				data[i] = static_cast<std::string>(token);
				myline.erase(0, pos + delimiter.length());
				++i;
			}
			if (nLine < nRecords) {
				if ('B' == data[1][0]) {
					target[nLine] = 0.0;
				}
				else {
					target[nLine] = 1.0;
				}
				for (int j = 0; j < nFeatures; ++j) {
					inputs[nLine][j] = atof(data[j + 2].c_str());
				}
				++nLine;
			}
		}

		myfile.close();
		return true;
	}
	else {
		return false;
	}
}

bool DataHolder::ReadDataSpam() {
	std::ifstream myfile;
	myfile.open("..//DATA//spambase//spambase.data");
	std::string delimiter = ",";
	nRecords = 4601;
	nFeatures = 57;

	target = std::make_unique<double[]>(nRecords);
	inputs = std::make_unique<std::unique_ptr<double[]>[]>(nRecords);
	for (int k = 0; k < nRecords; ++k) {
		inputs[k] = std::make_unique<double[]>(nFeatures);
	}

	std::string myline;
	auto data = std::make_unique<std::string[]>(59);
	int nLine = 0;
	if (myfile.is_open()) {
		while (myfile.good()) {
			myfile >> myline;
			myline += ",";
			size_t pos = 0;
			std::string token;
			int i = 0;
			while ((pos = myline.find(delimiter)) != std::string::npos) {
				token = myline.substr(0, pos);
				data[i] = static_cast<std::string>(token);
				myline.erase(0, pos + delimiter.length());
				++i;
			}
			if (nLine < nRecords) {
				for (int j = 0; j < nFeatures; ++j) {
					inputs[nLine][j] = atof(data[j].c_str());
				}
				target[nLine] = atof(data[nFeatures].c_str());
				++nLine;
			}
		}
		myfile.close();
		return true;
	}
	else {
		return false;
	}
}

bool DataHolder::ReadDataMFormula() {
	nRecords = 10000;
	nFeatures = 5;
	double pi = 3.14159265359;
	target = std::make_unique<double[]>(nRecords);
	inputs = std::make_unique<std::unique_ptr<double[]>[]>(nRecords);

	int counter = 0;
	while (true) {
		inputs[counter] = std::make_unique<double[]>(nFeatures);

		for (int i = 0; i < nFeatures; ++i) {
			inputs[counter][i] = static_cast<double>((rand() % 1000) / 1000.0);
		}

		double y = (1.0 / pi);
		y *= (2.0 + 2.0 * inputs[counter][2]);
		y *= (1.0 / 3.0);
		y *= atan(20.0 * (inputs[counter][0] - 0.5 + inputs[counter][1] / 6.0) * exp(inputs[counter][4])) + pi / 2.0;

		double z = (1.0 / pi);
		z *= (2.0 + 2.0 * inputs[counter][3]);
		z *= (1.0 / 3.0);
		z *= atan(20.0 * (inputs[counter][0] - 0.5 - inputs[counter][1] / 6.0) * exp(inputs[counter][4])) + pi / 2.0;

		target[counter] = static_cast<double>(y + z);
		if (++counter >= nRecords) break;
	}

	return true;
}

