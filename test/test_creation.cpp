#include <iostream>
#include <vt_linalg>

using namespace vt;

#define NEWLINE() std::cout<<'\n'

int main() {
    constexpr size_t N = 128;
    auto M1 = vt::make_numeric_matrix(
            {{1, 2, 3},
             {4, 5, 6}}
    );

    for (auto &row: M1) {
        for (auto &x: row) {
            std::cout << x << ' ';
        }
        std::cout << '\n';
    }

    auto vs = vt::detail::make_nested(
            {{1, 2, 3},
             {4, 5, 6}}
    );

    for (auto &row: vs) {
        for (auto &x: row) {
            std::cout << x << ' ';
        }
        std::cout << '\n';
    }

    return 0;
}
