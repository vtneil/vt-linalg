#include <iostream>
#include "standard_utility.h"
#include <vt_linalg>


using namespace vt;


int main() {
    numeric_matrix_static_t<numeric_matrix<2>, 2, 2> A;
    numeric_matrix_static_t<numeric_matrix<2>, 2, 2> B;
    A[0][0] = numeric_matrix<2, 2>::identity();
//    A[0][1] = numeric_matrix<2, 2>::identity();
//    A[1][0] = numeric_matrix<2, 2>::identity();
    A[1][1] = numeric_matrix<2, 2>::identity();

    B[0][0] = make_numeric_matrix({{1, 2},
                                   {3, 4}});
    B[0][1] = make_numeric_matrix({{1, 2},
                                   {3, 4}});
    B[1][0] = make_numeric_matrix({{1, 2},
                                   {3, 4}});
    B[1][1] = make_numeric_matrix({{1, 2},
                                   {3, 4}});

    auto C = A * B;
    auto unrolled_C = make_quad_matrix(C[0][0], C[0][1], C[1][0], C[1][1]);

    for (auto &row: unrolled_C) {
        for (auto &x: row) {
            std::cout << x << ' ';
        }
        std::cout << '\n';
    }

    auto v = make_numeric_vector(
            make_numeric_vector({1, 2, 3}),
            make_numeric_vector({4, 5, 6, 7}),
            make_numeric_vector({8, 9, 10}),
            make_numeric_vector({11, 12, 13})
    );

    for (auto &x: v) {
        std::cout << x << ' ';
    }
    std::cout << '\n';

    return 0;
}