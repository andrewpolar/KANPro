//Concept: Andrew Polar and Mike Poluektov
//Developer Andrew Polar

//License
//In case if end user finds the way of making a profit by using this code and earns
//billions of US dollars and meet developer bagging change in the street near McDonalds,
//he or she is not in obligation to buy him a sandwich.

//Symmetricity
//In case developer became rich and famous by publishing this code and meet misfortunate
//end user who went bankrupt by using this code, he is also not in obligation to buy
//end user a sandwich.

//Publications:
//https://www.sciencedirect.com/science/article/abs/pii/S0016003220301149
//https://www.sciencedirect.com/science/article/abs/pii/S0952197620303742
//https://arxiv.org/abs/2305.08194

#include <iostream>
#include "DataHolder.h"
#include "KANAddendPL.h"
#include "Basis.h"
#include "KANAddend.h"

void ShowMatrix(std::unique_ptr<std::unique_ptr<double[]>[]>& matrix, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            printf("%5.3f ", matrix[i][j]);
        }
        printf("\n");
    }
}

void ShowVector(std::unique_ptr<double[]>& ptr, int N) {
    int cnt = 0;
    for (int i = 0; i < N; ++i) {
        printf("%5.2f ", ptr[i]);
        if (++cnt >= 10) {
            printf("\n");
            cnt = 0;
        }
    }
}

void SwapRows(std::unique_ptr<double[]>& row1, std::unique_ptr<double[]>& row2, int cols) {
    auto ptr = std::make_unique<double[]>(cols);
    for (int i = 0; i < cols; ++i) {
        ptr[i] = row1[i];
    }
    for (int i = 0; i < cols; ++i) {
        row1[i] = row2[i];
    }
    for (int i = 0; i < cols; ++i) {
        row2[i] = ptr[i];
    }
}

void SwapScalars(double& x1, double& x2) {
    double buff = x1;
    x1 = x2;
    x2 = buff;
}

void Shuffle(std::unique_ptr<std::unique_ptr<double[]>[]>& matrix, std::unique_ptr<double[]>& vector, int rows, int cols) {
    for (int i = 0; i < 2 * rows; ++i) {
        int n1 = rand() % rows;
        int n2 = rand() % rows;
        SwapRows(matrix[n1], matrix[n2], cols);
        SwapScalars(vector[n1], vector[n2]);
    }
}

void FindMinMax(std::unique_ptr<double[]>& xmin, std::unique_ptr<double[]>& xmax,
    double& targetMin, double& targetMax,
    std::unique_ptr<std::unique_ptr<double[]>[]>& matrix,
    std::unique_ptr<double[]>& target, int nRows, int nCols) {

    for (int i = 0; i < nCols; ++i) {
        xmin[i] = DBL_MAX;
        xmax[i] = -DBL_MIN;
    }

    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            if (matrix[i][j] < xmin[j]) xmin[j] = static_cast<double>(matrix[i][j]);
            if (matrix[i][j] > xmax[j]) xmax[j] = static_cast<double>(matrix[i][j]);
        }
    }

    targetMin = DBL_MAX;
    targetMax = -DBL_MIN;
    for (int j = 0; j < nRows; ++j) {
        if (target[j] < targetMin) targetMin = target[j];
        if (target[j] > targetMax) targetMax = target[j];
    }
}

int main() {

    srand((unsigned int)time(NULL));
    auto dataHolder = std::make_unique<DataHolder>();
    bool status = dataHolder->ReadDataMFormula();
    if (false == status) {
        printf("Failed to open file");
        exit(0);
    }

    clock_t start_application = clock();
    Shuffle(dataHolder->inputs, dataHolder->target, dataHolder->nRecords, dataHolder->nFeatures);
    
    //ShowMatrix(dataHolder->inputs, dataHolder->nRecords, dataHolder->nFeatures);
    //ShowVector(dataHolder->target, dataHolder->nRecords);

    auto xmin = std::make_unique<double[]>(dataHolder->nFeatures);
    auto xmax = std::make_unique<double[]>(dataHolder->nFeatures);
    double targetMin;
    double targetMax;

    FindMinMax(xmin, xmax, targetMin, targetMax, dataHolder->inputs, dataHolder->target, 
        dataHolder->nRecords, dataHolder->nFeatures);

    //Initialization of KAN and training for Formula3
    int nRecords = dataHolder->nRecords;   //10 000
    int nFeatures = dataHolder->nFeatures; //5
    int nModels = 11;
    int PWLEpochs = 40;
    int SplineEpochs = 50;
    int SplineTopComplete = 6;
    int SplineBtmComplete = 6;
    double sensitivity = 0.01 * (targetMax - targetMin);
    int innerPoints = 2;
    int outerPoints = 2;
    double muInnerPL = 0.01;
    double muOuterPL = 0.01;
    double muInnerSP = 0.05;
    double muOuterSP = 0.05;

    //initialization of piecewise linear model
    auto addends = std::make_unique<std::unique_ptr<KANAddendPL>[]>(nModels);
    for (int i = 0; i < nModels; ++i) {
        addends[i] = std::make_unique<KANAddendPL>(xmin, xmax, targetMin / nModels, targetMax / nModels, innerPoints,
            outerPoints, muInnerPL, muOuterPL, nFeatures);
    }

    //split data into training and validation
    auto isTraining = std::make_unique<bool[]>(nRecords);
    for (int i = 0; i < nRecords; ++i) {
        if (rand() % 1000 < 620) {
            isTraining[i] = true;
        }
    }
 
    //training of piecewise linear model with incrementing of linear segments
    for (int epoch = 0; epoch < PWLEpochs; ++epoch) {
        for (int i = 0; i < nRecords; ++i) {
            if (!isTraining[i]) continue;
            double model = 0.0;
            for (int j = 0; j < nModels; ++j) {
                model += addends[j]->ComputeUsingInput(dataHolder->inputs[i]);
            }
            double residual = dataHolder->target[i] - model;
            for (int j = 0; j < nModels; ++j) {
                addends[j]->UpdateUsingMemory(residual);
            }
        }
        printf("PWL training epoch %d\r", epoch);

        if (epoch > 0 && epoch < 9) {
            for (int j = 0; j < nModels; ++j) {
                addends[j]->IncrementInner();
            }
        }

        if (epoch > 0 && epoch < 13) {
            for (int j = 0; j < nModels; ++j) {
                addends[j]->IncrementOuter();
            }
        }
    }
    printf("\n");
    clock_t end_PWL_training = clock();
    printf("Time for PWL training %2.3f sec.\n", (double)(end_PWL_training - start_application) / CLOCKS_PER_SEC);

    //initialization of spline model
    int innerP = addends[0]->HowManyInner();
    int outerP = addends[0]->HowManyOuter();
 
    //Next part is converting piecewise linear model into splines
    Basis innerBasis(innerP);
    Basis outerBasis(outerP);

    std::vector<std::unique_ptr<KANAddend>> splineaddends;
    for (int i = 0; i < nModels; ++i) {
        splineaddends.push_back(std::make_unique<KANAddend>(xmin, xmax,
            targetMin / nModels, targetMax / nModels, innerBasis, outerBasis, muInnerSP, muOuterSP, nFeatures));
    }

    //making spline model as a copy of piecewise linear
    for (int i = 0; i < nModels; ++i) {
        auto x = addends[i]->GetAllOuterPoints();
        splineaddends[i]->UpdateOuterPoints(x, outerP);
        for (int j = 0; j < nFeatures; ++j) {
            auto x = addends[i]->_u->GetUPoints(j);
            splineaddends[i]->_u->AssignUPoints(x, j, innerP);
        }
    }

    //next block is training of spline model
    auto residuals = std::make_unique<double[]>(nRecords);
    for (int step = 0; step < SplineEpochs; ++step) {
        for (int i = 0; i < nRecords; ++i) {
            if (!isTraining[i]) continue;
            if (step > SplineTopComplete && step < SplineEpochs - SplineBtmComplete && residuals[i] < sensitivity) continue;
            double model = 0.0;
            for (int j = 0; j < nModels; ++j) {
                model += splineaddends[j]->ComputeUsingInput(dataHolder->inputs[i]);
            }
            double diff = dataHolder->target[i] - model;
            for (int j = 0; j < nModels; ++j) {
                splineaddends[j]->UpdateUsingMemory(diff);
            }
            residuals[i] = fabs(diff);
        }
        printf("Spline training epoch %d\r", step);
    }
    printf("\n");

    //computing residual error for validation sample
    double error = 0.0;
    int cnt = 0;
    for (int i = 0; i < nRecords; ++i) {
        if (isTraining[i]) continue;
        double model = 0.0;
        for (int j = 0; j < nModels; ++j) {
            model += splineaddends[j]->ComputeUsingInput(dataHolder->inputs[i]);
        }
        double diff = dataHolder->target[i] - model;
        error += diff * diff;
        ++cnt;
    }
    error /= cnt;
    error = sqrt(error);
    error /= targetMax - targetMin;

    printf("The realive RMSE for validation sample %6.4f\n", error);

    clock_t end_encoding = clock();
    printf("Time for training %2.3f sec.\n\n", (double)(end_encoding - start_application) / CLOCKS_PER_SEC);
}
