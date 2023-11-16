#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

class Matrix {
private:
    std::vector<std::vector<double>> data;
public:
    Matrix(const std::vector<std::vector<double>>& data);
    bool EqMatrix(const Matrix& other) const;
    void SumMatrix(const Matrix& other);
    void SubMatrix(const Matrix& other);
    void MulNumber(double number);
    void MulMatrix(const Matrix& other);
    double Determinant() const;
    Matrix Transpose() const;
    Matrix CalcComplements() const;
    Matrix InverseMatrix() const;
};

#endif
