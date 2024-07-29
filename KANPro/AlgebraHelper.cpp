#include "AlgebraHelper.h"
#include <vector>

using namespace std;

void AlgebraHelper::ShowMatrix(std::unique_ptr<std::unique_ptr<double[]>[]>& matrix, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            printf("%5.3f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void AlgebraHelper::MultiplySquare(std::unique_ptr<std::unique_ptr<double[]>[]>& X, std::unique_ptr<std::unique_ptr<double[]>[]>& Y,
    std::unique_ptr<std::unique_ptr<double[]>[]>& Z, int size) {

    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            Z[i][j] = 0.0;
            for (int k = 0; k < size; ++k)
            {
                Z[i][j] += X[i][k] * Y[k][j];
            }
        }
    }
}

bool AlgebraHelper::choldc1(int n, std::unique_ptr<std::unique_ptr<double[]>[]>& a, std::unique_ptr<double[]>& p) {
    int i, j, k;
    double sum;

    for (i = 0; i < n; i++)
    {
        for (j = i; j < n; j++)
        {
            sum = a[i][j];
            for (k = i - 1; k >= 0; k--)
            {
                sum -= a[i][k] * a[j][k];
            }
            if (i == j)
            {
                if (sum <= 0)
                {
                    printf("%s\n", "Matrix is not positive definite");
                    return false;
                }
                p[i] = sqrt(sum);
            }
            else
            {
                a[j][i] = sum / p[i];
            }
        }
    }
    return true;
}

bool AlgebraHelper::choldc(int n, std::unique_ptr<std::unique_ptr<double[]>[]>& A, std::unique_ptr<std::unique_ptr<double[]>[]>& a) {
    int i, j;
    auto p = std::make_unique<double[]>(n);

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            a[i][j] = A[i][j];

    if (!choldc1(n, a, p))
    {
        return false;
    }

    for (i = 0; i < n; i++)
    {
        a[i][i] = p[i];
        for (j = i + 1; j < n; j++)
        {
            a[i][j] = 0;
        }
    }

    return true;
}

bool AlgebraHelper::CholeskySolution(std::unique_ptr<std::unique_ptr<double[]>[]>& M, std::unique_ptr<std::unique_ptr<double[]>[]>& Inv, int size) {
    auto T = std::make_unique<std::unique_ptr<double[]>[]>(size);
    for (int k = 0; k < size; ++k) {
        T[k] = std::make_unique<double[]>(size);
    }

    //find lower triangular
    if (!choldc(size, M, T))
    {
        return false;
    }
    //make it symmetric
    for (int i = 0; i < size; ++i)
        for (int j = i + 1; j < size; ++j)
            T[i][j] = T[j][i];

    auto I = std::make_unique<std::unique_ptr<double[]>[]>(size);
    for (int k = 0; k < size; ++k) {
        I[k] = std::make_unique<double[]>(size);
    }
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            if (i == j) I[i][j] = 1.0;
            else I[i][j] = 0.0;
        }
    }
    auto v = std::make_unique<double[]>(size);

    for (int k = 0; k < size; ++k)
    {
        for (int i = 0; i < size; ++i)
        {
            v[i] = 0;
        }
        //solve lower triangular
        for (int i = 0; i < size; ++i)
        {
            double sum = 0.0;
            for (int j = 0; j < i; ++j)
            {
                sum += T[i][j] * v[j];
            }
            v[i] = (I[k][i] - sum) / T[i][i];
        }
        //solve upper triangular
        for (int i = size - 1; i >= 0; --i)
        {
            double sum = 0.0;
            for (int j = size - 1; j > i; --j)
            {
                sum += T[i][j] * Inv[k][j];
            }
            Inv[k][i] = (v[i] - sum) / T[i][i];
        }
    }
    return true;
}

void AlgebraHelper::InverseSelfTest() {
    const int size = 4;
    auto H = std::make_unique<std::unique_ptr<double[]>[]>(size);
    for (int k = 0; k < size; ++k) {
        H[k] = std::make_unique<double[]>(size);
    }
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            H[i][j] = 1.0 / (double)(i + j + 1);
        }
    }

    auto Inv = std::make_unique<std::unique_ptr<double[]>[]>(size);
    for (int k = 0; k < size; ++k) {
        Inv[k] = std::make_unique<double[]>(size);
    }

    if (CholeskySolution(H, Inv, size))
    {
        ShowMatrix(Inv, size, size);
        auto I = std::make_unique<std::unique_ptr<double[]>[]>(size);
        for (int k = 0; k < size; ++k) {
            I[k] = std::make_unique<double[]>(size);
        }
        MultiplySquare(Inv, H, I, size);
        ShowMatrix(I, size, size);
    }
    else
    {
        printf("Failed to obtain cholessky solution\n");
    }

    if (gaussJordan(H, Inv, size))
    {
        ShowMatrix(Inv, size, size);
        auto I = std::make_unique<std::unique_ptr<double[]>[]>(size);
        for (int k = 0; k < size; ++k) {
            I[k] = std::make_unique<double[]>(size);
        }
        MultiplySquare(Inv, H, I, size);
        ShowMatrix(I, size, size);
    }
    else
    {
        printf("Failed to obtain GaussJordan solution\n");
    }
}

bool AlgebraHelper::gaussJordan(std::unique_ptr<std::unique_ptr<double[]>[]>& matrix, std::unique_ptr<std::unique_ptr<double[]>[]>& inverse, int size) {
    int n = size;
    vector<vector<double>> augmentedMatrix(n, vector<double>(2 * n));

    // Create the augmented matrix [matrix | I]
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            augmentedMatrix[i][j] = matrix[i][j];
            augmentedMatrix[i][j + n] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Perform row operations to convert [matrix | I] to [I | inverse]
    for (int i = 0; i < n; ++i) {
        // Find pivot
        double pivot = augmentedMatrix[i][i];
        if (pivot == 0) {
            return false; // Matrix is singular
        }

        // Normalize the pivot row
        for (int j = 0; j < 2 * n; ++j) {
            augmentedMatrix[i][j] /= pivot;
        }

        // Eliminate the column
        for (int k = 0; k < n; ++k) {
            if (k != i) {
                double factor = augmentedMatrix[k][i];
                for (int j = 0; j < 2 * n; ++j) {
                    augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                }
            }
        }
    }

    // Extract the inverse matrix
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            inverse[i][j] = augmentedMatrix[i][j + n];
        }
    }

    return true;
}

//https://people.clas.ufl.edu/kees/files/CubicSplines.pdf
std::unique_ptr<std::unique_ptr<double[]>[]> AlgebraHelper::GenerateTriDiagonal(int N, std::unique_ptr<double[]>& h) {
    auto M = std::make_unique<std::unique_ptr<double[]>[]>(N);
    for (int k = 0; k < N; ++k) {
        M[k] = std::make_unique<double[]>(N);
    }

    M[0][0] = 1.0;
    for (int j = 1; j < N; ++j)
    {
        M[0][j] = 0.0;
    }
    for (int i = 1; i < N - 1; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if (i == j) M[i][j] = 2.0 * (h[i - 1] + h[i]);
            else if (1 == i - j) M[i][j] = h[i - 1];
            else if (1 == j - i) M[i][j] = h[j - 1];
            else M[i][j] = 0.0;
        }
    }
    for (int j = 0; j < N - 1; ++j)
    {
        M[N - 1][j] = 0.0;
    }
    M[N - 1][N - 1] = 1.0;

    return M;
}

//https://people.clas.ufl.edu/kees/files/CubicSplines.pdf
void AlgebraHelper::MakeSplines(std::unique_ptr<std::unique_ptr<double[]>[]>& A,
    std::unique_ptr<double[]>& y, std::unique_ptr<double[]>& h,
    std::unique_ptr<double[]>& a, std::unique_ptr<double[]>& b, std::unique_ptr<double[]>& c, std::unique_ptr<double[]>& d, int lenY)
{
    int N = lenY;

    auto z = std::make_unique<double[]>(N);
    z[0] = 0.0;
    for (int i = 1; i < N - 1; ++i)
    {
        z[i] = 3.0 * (y[i + 1] - y[i]) / h[i] - 3.0 * (y[i] - y[i - 1]) / h[i - 1];
    }
    z[N - 1] = 0.0;

    auto v = std::make_unique<double[]>(N);
    for (int i = 0; i < N; ++i)
    {
        v[i] = 0.0;
        for (int j = 0; j < N; ++j)
        {
            v[i] += A[i][j] * z[j];
        }
    }

    //auto a = std::make_unique<double[]>(N - 1);
    for (int i = 0; i < N - 1; ++i)
    {
        a[i] = y[i];
    }

    //auto b = std::make_unique<double[]>(N - 1);
    for (int i = 0; i < N - 1; ++i)
    {
        b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (2.0 * v[i] + v[i + 1]) / 3.0;
    }

    //auto c = std::make_unique<double[]>(N - 1);
    for (int i = 0; i < N - 1; ++i)
    {
        c[i] = v[i];
    }

    //auto d = std::make_unique<double[]>(N - 1);
    for (int i = 0; i < N - 1; ++i)
    {
        d[i] = (v[i + 1] - v[i]) / 3.0 / h[i];
    }
}

void AlgebraHelper::SplinesSelfTest()
{
    const int lenXY = 4;
    const int lenh = 3;
    double x[lenXY] = { 0.0, 1.0, 2.0, 2.5 };
    auto y = std::make_unique<double[]>(lenXY);
    y[0] = 0.0;
    y[1] = 1.0;
    y[2] = 8.0;
    y[3] = 9.0;

    auto h = std::make_unique<double[]>(lenh);
    for (int i = 0; i < lenh; ++i)
    {
        h[i] = x[i + 1] - x[i];
    }
    auto M = GenerateTriDiagonal(lenXY, h);

    auto R = std::make_unique<std::unique_ptr<double[]>[]>(lenXY);
    for (int k = 0; k < lenXY; ++k) {
        R[k] = std::make_unique<double[]>(lenXY);
    }

    gaussJordan(M, R, lenXY);

    auto a = std::make_unique<double[]>(lenXY);
    auto b = std::make_unique<double[]>(lenXY);
    auto c = std::make_unique<double[]>(lenXY);
    auto d = std::make_unique<double[]>(lenXY);

    MakeSplines(R, y, h, a, b, c, d, lenXY);

    for (int j = 0; j < lenh; ++j)
    {
        printf("%4.2f %4.2f %4.2f %4.2f\n", a[j], b[j], c[j], d[j]);
    }
    printf("\n");

    //expected result
    //0.000 - 1.091 -0.000  2.091
    //1.000   5.182  6.273 -4.455
    //8.000   4.364 -7.091  4.727
}

