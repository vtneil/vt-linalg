#include "vnet_linalg/vector.h"
#include "vnet_linalg/matrix.h"
#include <iostream>

int main() {
    Vector<int> v1 = Vector<int>::from(1, 2, 3, 4);
    Vector<int> v2 = Vector<int>::from(5, 6, 7, 1);

    Matrix<int> m1(2, 3, 6);
    Matrix<int> m2(m1);

    m2[0][0] = 2;
    m2[1][1] = 3;

    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            std::cout << m2[i][j] << ' ';
        }
        std::cout << '\n';
    }

    for (auto &x: m2.diag()) {
        std::cout << x << ' ';
    }

    return 0;
}
