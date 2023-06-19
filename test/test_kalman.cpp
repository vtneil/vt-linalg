#include <iostream>
#include "kalman.h"

using namespace vt;

real_t dt = 0.5;
real_t base_noise_value = 0.001;
numeric_matrix<2, 2> F({{1, dt},
                        {0, 1}});
numeric_matrix<2, 1> B({{0},
                        {0}});
numeric_matrix<1, 2> H({{1, 0}});
numeric_matrix<2, 2> Q = numeric_matrix<2, 2>::diagonals(base_noise_value);
numeric_matrix<1, 1> R = numeric_matrix<1, 1>::diagonals(base_noise_value);
numeric_vector<2> x0({0, 0});
kalman_filter_t<2, 1, 1> kf(F, B, H, Q, R, x0);

int main() {
    constexpr size_t SAMPLES = 16;
    numeric_vector<1> controls[SAMPLES] = {};
    numeric_vector<1> measurements[SAMPLES] = {};

    real_t i_r = 0;
    for (auto &x: measurements) x = {i_r++};

    for (size_t i = 0; i < SAMPLES; ++i) {
        kf.predict(controls[i]);
        std::cout << "Predicted state: ";
        for (const auto &state: kf.state_vector()) {
            std::cout << state << " ";
        }
        std::cout << '\n';

        kf.update(measurements[i]);
        std::cout << i + 1 << " - Estimated state: ";
        for (const auto &state: kf.state_vector()) {
            std::cout << state << " ";
        }
        std::cout << '\n';
    }

    return 0;
}
