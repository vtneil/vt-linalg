#include <iostream>
#include "../vt_linalg"

#define NEWLINE() std::cout<<'\n'


int main() {
    numeric_matrix<1, 3> A({{3, 4, 2}});
    numeric_matrix<3, 4> B({{13, 9, 7, 15},
                            {8,  7, 4, 6},
                            {6,  4, 0, 3}});

    auto C = A * B;
    for (auto &row: C) {
        for (auto &x: row) {
            std::cout << x << ' ';
        }
        std::cout << '\n';
    }
    return 0;
}
