#include <iostream>
#include <vector>

// Функция для генерации исключения при неправильном размере матрицы
void checkMatrixSize(const std::vector<std::vector<int>>& matrix1, const std::vector<std::vector<int>>& matrix2) {
    if (matrix1.size() != matrix2.size() || matrix1[0].size() != matrix2[0].size()) {
        throw "Матрицы должны быть одинакового размера";
    }
}

// Функция для сложения матриц
std::vector<std::vector<int>> addMatrix(const std::vector<std::vector<int>>& matrix1, const std::vector<std::vector<int>>& matrix2) {
    checkMatrixSize(matrix1, matrix2); // Проверка размера матриц
    std::vector<std::vector<int>> result(matrix1.size(), std::vector<int>(matrix1[0].size()));

    for (int i = 0; i < matrix1.size(); i++) {
        for (int j = 0; j < matrix1[0].size(); j++) {
            result[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
    return result;
}

// Функция для вычитания матриц
std::vector<std::vector<int>> subtractMatrix(const std::vector<std::vector<int>>& matrix1, const std::vector<std::vector<int>>& matrix2) {
    checkMatrixSize(matrix1, matrix2); // Проверка размера матриц
    std::vector<std::vector<int>> result(matrix1.size(), std::vector<int>(matrix1[0].size()));

    for (int i = 0; i < matrix1.size(); i++) {
        for (int j = 0; j < matrix1[0].size(); j++) {
            result[i][j] = matrix1[i][j] - matrix2[i][j];
        }
    }
    return result;
}

// Функция для умножения матрицы на скаляр
std::vector<std::vector<int>> multiplyMatrixByScalar(const std::vector<std::vector<int>>& matrix, const int scalar) {
    std::vector<std::vector<int>> result(matrix.size(), std::vector<int>(matrix[0].size()));

    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[0].size(); j++) {
            result[i][j] = matrix[i][j] * scalar;
        }
    }
    return result;
}

// Функция для перемножения матриц
std::vector<std::vector<int>> multiplyMatrix(const std::vector<std::vector<int>>& matrix1, const std::vector<std::vector<int>>& matrix2) {
    if (matrix1[0].size() != matrix2.size()) {
        throw "Количество столбцов первой матрицы должно быть равно количеству строк второй матрицы";
    }

    std::vector<std::vector<int>> result(matrix1.size(), std::vector<int>(matrix2[0].size()));

    for (int i = 0; i < matrix1.size(); i++) {
        for (int j = 0; j < matrix2[0].size(); j++) {
            result[i][j] = 0;
            for (int k = 0; k < matrix1[0].size(); k++) {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
    return result;
}

// Функция для транспонирования матрицы
std::vector<std::vector<int>> transposeMatrix(const std::vector<std::vector<int>>& matrix) {
    std::vector<std::vector<int>> result(matrix[0].size(), std::vector<int>(matrix.size()));

    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[0].size(); j++) {
            result[j][i] = matrix[i][j];
        }
    }
    return result;
}
