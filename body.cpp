#include <iostream>
#include "matrix.h"

Matrix::Matrix(const std::vector<std::vector<double>>& data) : data(data) {}

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

void Matrix::Fill {
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j = cols_; j++) {
            matrix_[i][j] = 0;
        }
    }
}

void Matrix::OutputMatrix() const {
    for (const auto& row : data) {
        for (const auto& elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
}

void Matrix::SumMatrix(const Matrix& other) {
    if (data.size() != other.data.size() || data[0].size() != other.data[0].size())
        throw std::invalid_argument("Matrices have different sizes!");

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

double** Matrix::AllocateMatrix(int rows, int cols) {
    double** matrix = new double[rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new double[cols];
    }
    return matrix;
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
    if (data.size() != data[0].size()) {
        throw std::invalid_argument("Matrix is not square!");
    }

    int n = data.size();
    std::vector<std::vector<double>> tempMatrix = data;
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
    int rows = data.size();
    int cols = data[0].size();

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
    for (int i = 0; i < transposed.data.size(); i++) {
        for (int j = 0; j < transposed.data[i].size(); j++) {
            transposed.data[i][j] /= det;
        }
    }

    return transposed;
}

void Print() const {
    for (const auto& row : data) {
        for (const auto& elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
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
