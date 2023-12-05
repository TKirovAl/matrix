#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

class Matrix {
private:
    int data, cols_, rows_;
    double** matrix_;
    double** AllocateMatrix(int rows, int cols);
    void Fill();
    
public:
    Matrix();
    Matrix(int rows, int cols);
    int getRows();
    int getCols();
    void OutputMatrix() noexcept;
    void Fill(double** matrix_);
    bool EqMatrix(const Matrix& other) const;
    void OutputMatrix() const;
    void SubMatrix(const Matrix& other);
    void SumMatrix(const Matrix& other);
    void MulMatrix(const Matrix& other);
    Matrix Transpose() const;
    void MulNumber(double number);
    double Determinant() const;
};

#endif
