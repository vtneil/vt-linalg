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
    Vector<double> v1 = Vector<double>({7, 2, 3, 1});
    Vector<double> v2 = Vector<double>({1, 2, 3, 4});

    Matrix<double> A = Matrix<double>({{1, 2, 3},
                                       {4, 5, 6},
                                       {7, 8, 10}});

    print_matrix(A);
    std::cout << "-----\n";

    print_matrix(RRE(A));
    std::cout << "-----\n";

    print_matrix(inv(A));
    std::cout << "-----\n";

    std::cout << tr(A) << '\n';
    std::cout << "-----\n";

    std::cout << det(A) << '\n';
    std::cout << "-----\n";

    return 0;
}
