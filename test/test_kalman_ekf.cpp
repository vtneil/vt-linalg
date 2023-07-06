#include <iostream>
#include "vt_kalman"

using namespace vt;

real_t dt = 1.0;
real_t hdts = 0.5 * dt * dt;
real_t base_noise_value = 0.001;

numeric_vector<3> f(const numeric_vector<3> &x,
                    const numeric_vector<1> &) {
    return numeric_vector<3>({x[0] + dt * x[1] + hdts * x[2],
                              x[1] + dt * x[2],
                              x[2]});
}

numeric_matrix<3, 3> Fj(const numeric_vector<3> &,
                        const numeric_vector<1> &) {
    return numeric_matrix<3, 3>({{1, dt, hdts},
                                 {0, 1,  dt},
                                 {0, 0,  1}});
}

numeric_vector<1> h(const numeric_vector<3> &x) {
    numeric_vector<1> z;
    z[0] = x[0];
    return z;
}

numeric_matrix<1, 3> Hj(const numeric_vector<3> &) {
    return numeric_matrix<1, 3>({{1}});
}

numeric_matrix<3, 3> Q = numeric_matrix<3, 3>::diagonals(base_noise_value);
numeric_matrix<1, 1> R = numeric_matrix<1, 1>::diagonals(base_noise_value);
numeric_vector<3> x0; // {x, v, a}

extended_kalman_filter_t<3, 1, 1> ekf(f, Fj, h, Hj, Q, R, x0);

int main() {
    constexpr size_t SAMPLES = 16;
    numeric_vector<1> controls[SAMPLES] = {};
    numeric_vector<1> measurements[SAMPLES] = {};

    real_t i_r = 0;
    for (auto &x: measurements) x = {i_r++};

    for (size_t i = 0; i < SAMPLES; ++i) {
        ekf.predict(controls[i]);
        std::cout << "Predicted state: ";
        for (const auto &state: ekf.state_vector) {
            std::cout << state << " ";
        }
        std::cout << '\n';

        ekf.update(measurements[i]);
        std::cout << i + 1 << " - Estimated state: ";
        for (const auto &state: ekf.state_vector) {
            std::cout << state << " ";
        }
        std::cout << '\n';
    }

    return 0;
}
