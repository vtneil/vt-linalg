#include <iostream>
#include "standard_utility.h"
#include <vt_linalg>


using namespace vt;


int main() {
    auto v1 = vt::detail::make_numeric_vector_static_t<real_t>({1, 2, 3}, {4, 5});
    numeric_matrix<2, 3> zero;
    numeric_matrix<2, 3> A1 = make_numeric_matrix({{1, 2, 3},
                                                   {4, 5, 6}});
    numeric_matrix<2, 3> A2 = make_numeric_matrix({{7,  8,  9},
                                                   {10, 11, 12}});
    numeric_matrix<2, 3> A3 = make_numeric_matrix({{3,  5,  7},
                                                   {-1, -2, -3}});
//    numeric_matrix<2, 3> A = numeric_matrix<2, 3>({{A1}});
    auto A = make_block_matrix(
            {{A1, A2, A3, A1},
             {A1, A2, A3, A1},
             {A1, A2, A3, A1}}
    );

    auto B = make_block_matrix({{A, A}}).transpose();

    for (auto &x: v1) {
        std::cout << x << ' ';
    }
    std::cout << '\n';

    for (auto &row: B) {
        for (auto &x: row) {
            std::cout << x << ' ';
        }
        std::cout << '\n';
    }

    return 0;
}