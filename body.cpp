#include <iostream>
#include "matrix.h"

Matrix::Matrix(){
    
}

bool Matrix::EqMatrix(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_)
        return false;
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            if (matrix_[i][j] != other.matrix_[i][j])
                return false;
        }
    }
    return true;
}

void Matrix::Fill() {
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j = cols_; j++) {
            matrix_[i][j] = 0;
        }
    }
}

void Matrix::OutputMatrix() const {
    
}

void Matrix::SumMatrix(const Matrix& other) {
    if (rows_ != other.rows_ || cols_ != other.cols_)
        throw std::invalid_argument("Matrices have different sizes!");

    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            matrix_[i][j] += other.matrix_[i][j];
        }
    }
}

void Matrix::SubMatrix(const Matrix& other) {
    if (rows_ != other.rows_ || cols_ != other.cols_)
        throw std::runtime_error("Matrices have different sizes!");

    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < rows_; j++) {
            matrix_[i][j] -= other.matrix_[i][j];
        }
    }
}

void Matrix::MulNumber(double number) {
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            matrix_[i][j] *= number;
        }
    }
}

double** Matrix::AllocateMatrix(int rows, int cols) {
    double** matrix = new double*[rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new double[cols];
    }
    return matrix;
}

void Matrix::MulMatrix(const Matrix& other) {
    if (cols_ != other.rows_)
        throw std::runtime_error("Incompatible matrix sizes for multiplication!");

    std::vector<std::vector<double>> result(rows_, std::vector<double>(other.cols_, 0));

    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < other.cols_; j++) {
            for (int k = 0; k < cols_; k++) {
                result[i][j] += matrix_[i][k] * other.matrix_[k][j];
            }
        }
    }

    matrix_ = result;
}

double Matrix::Determinant() const {
    if (rows_ != cols_) {
        throw std::invalid_argument("Matrix is not square!");
    }

    int n = rows_;
    std::vector<std::vector<double>> tempMatrix = matrix_;
    double det = 1;

    for (int i = 0; i < n; i++) {
        if (tempMatrix[i][i] == 0) {
            bool swapSuccess = false;
            for (int j = i + 1; j < n; j++) {
                if (tempMatrix[j][i] != 0) {
                    std::swap(tempMatrix[i], tempMatrix[j]);
                    det *= -1;  // меняем знак определителя при перестановке строк
                    swapSuccess = true;
                    break;
                }
            }
            if (!swapSuccess) {
                return 0;  // если все элементы в столбце равны нулю, определитель равен 0
            }
        }

        for (int j = i + 1; j < n; j++) {
            double coef = tempMatrix[j][i] / tempMatrix[i][i];
            for (int k = i; k < n; k++) {
                tempMatrix[j][k] -= tempMatrix[i][k] * coef;
            }
        }
    }
    for (int i = 0; i < n; i++) {
        det *= tempMatrix[i][i];
    }
    return det;
}

Matrix::Transpose() const {
    int rows = rows_;
    int cols = cols_);

    std::vector<std::vector<double>> transposed(cols, std::vector<double>(rows, 0));

    for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
            transposed[j][i] = data[i][j];
        }
    }

    return Matrix(transposed);
}

Matrix::InverseMatrix() const {
    double det = Determinant();
    if (det == 0) {
        throw std::runtime_error("The matrix is singular and does not have an inverse.");
    }

    Matrix::complements = CalcComplements();
    Matrix::transposed = complements.Transpose();
    for (int i = 0; i < transposed.rows_; i++) {
        for (int j = 0; j < transposed.cols_; j++) {
            transposed.matrix_[i][j] /= det;
        }
    }

    return transposed;
}

int main() {
    
    try {
        std::vector<std::vector<double>> matData = {{1, 2}, {3, 4}};
        Matrix::Matrix(matData);
        
        double det = mat.Determinant();
        Matrix transposed = mat.Transpose();
        Matrix complements = mat.CalcComplements();
        Matrix inverse = mat.InverseMatrix();
        Matrix subMatrix = mat.SubMatrix();
        Matrix sumMatrix = mat.SumMatrix();
        Matrix mulMatrix = mat.MulMatrix();
        
    }
    catch(const std::exception& e) {
        std::cout << "Ошибка: " << e.what() << std::endl;
    }
    
    return 0;
}
