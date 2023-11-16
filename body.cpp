#include <iostream>
#include "matrix.h"

Matrix::Matrix(const std::vector<std::vector<double>>& data) : data(data) {}

bool Matrix::EqMatrix(const Matrix& other) const {
    if (data.size() != other.data.size() || data[0].size() != other.data[0].size())
        return false;

    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[i].size(); j++) {
            if (data[i][j] != other.data[i][j])
                return false;
        }
    }

    return true;
}

void Matrix::SumMatrix(const Matrix& other) {
    if (data.size() != other.data.size() || data[0].size() != other.data[0].size())
        throw std::runtime_error("Matrices have different sizes!");

    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[i].size(); j++) {
            data[i][j] += other.data[i][j];
        }
    }
}

void Matrix::SubMatrix(const Matrix& other) {
    if (data.size() != other.data.size() || data[0].size() != other.data[0].size())
        throw std::runtime_error("Matrices have different sizes!");

    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[i].size(); j++) {
            data[i][j] -= other.data[i][j];
        }
    }
}

void Matrix::MulNumber(double number) {
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[i].size(); j++) {
            data[i][j] *= number;
        }
    }
}

void Matrix::MulMatrix(const Matrix& other) {
    if (data[0].size() != other.data.size())
        throw std::runtime_error("Incompatible matrix sizes for multiplication!");

    std::vector<std::vector<double>> result(data.size(), std::vector<double>(other.data[0].size(), 0));

    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < other.data[0].size(); j++) {
            for (int k = 0; k < data[0].size(); k++) {
                result[i][j] += data[i][k] * other.data[k][j];
            }
        }
    }

    data = result;
}

double Matrix::Determinant() const {
    if (data.size() != data[0].size())
        throw std::runtime_error("Matrix is not square!");

    int n = data.size();

    if (n == 1)
        return data[0][0];

    double determinant = 0;

    for (int i = 0; i < n; i++) {
        std::vector<std::vector<double>> submatrix(n - 1, std::vector<double>(n - 1));

        for (int j = 1; j < n; j++) {
            int k = 0;
            for (int l = 0; l < n; l++) {
                if (l == i)
                    continue;
                submatrix[j - 1][k] = data[j][l];
                k++;
            }
        }

      //  determinant += data[0][i] * (
