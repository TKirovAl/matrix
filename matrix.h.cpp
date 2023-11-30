#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

class Matrix {
private:
    int data, cols_;
    double** matrix_;
    double** AllocateMatrix(int rows, int cols);
    void Fill();
    
public:
    Matrix() noexcept;
    Matrix(int rows, int cols);
    int getRows();
    int getCols();
    void OutputMatrix() noexcept;
    void Fill(double** matrix_);
};

#endif
