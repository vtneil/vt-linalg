#include <iostream>
#include "../src/kalman.h"

using namespace vt;

numeric_matrix<2, 2> F({{1, 1},
                        {0, 1}});
numeric_matrix<1, 2> H({{1, 0}});
numeric_matrix<2, 2> Q({{0.001, 0},
                        {0,     0.001}});
numeric_matrix<1, 1> R({{0.01}});
numeric_vector<2> x0({0, 0});
simple_kalman_filter_t<2, 1> kf(F, H, Q, R, x0);

int main() {
    numeric_vector<1> measurements[] = {
            make_numeric_vector({1}),
            make_numeric_vector({2}),
            make_numeric_vector({3}),
            make_numeric_vector({4}),
            make_numeric_vector({5}),
    };

    for (const auto &measurement: measurements) {
        kf.predict();
        kf.update(measurement);

        std::cout << "Estimated state: ";
        for (const auto &state: kf.state_vector()) {
            std::cout << state << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
