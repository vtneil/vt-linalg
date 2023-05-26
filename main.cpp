#include "vnet_linalg/vector.h"
#include "vnet_linalg/matrix.h"
#include <iostream>

template<typename T>
void print_matrix(const Matrix<T> &A) {
    for (auto &x: A) {
        for (auto &y: *x) {
            std::cout << y << ' ';
        }
        std::cout << '\n';
    }
}

template<typename T>
void print_vector(const Vector<T> &v) {
    for (auto &y: v) {
        std::cout << y << ' ';
    }
    std::cout << '\n';
}

int main() {
    Vector<double> v1 = Vector<double>::from(7, 2, 3, 1);
    Vector<double> v2 = Vector<double>::from(3, 2, 1);

    Matrix<double> A = Matrix<double>::from(
            Vector<double>::from(4, 3),
            Vector<double>::from(6, 3)
    );

    Matrix<double> C = Matrix<double>::from_col(Vector<double>::from(1, 2, 3, 4));

    print_matrix((A ^ 2) + (2.0 * A) + Matrix<double>::identity(2));

    return 0;
}
