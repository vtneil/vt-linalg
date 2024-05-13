#include "vt_imu"
#include "vt_kalman"
#include <iostream>

using namespace vt;

constexpr real_t dt               = 1.0;
constexpr real_t hdts             = 0.5 * dt * dt;
constexpr real_t base_noise_value = 0.001;

constexpr numeric_matrix<3, 3> F({{1, dt, hdts},
                                  {0, 1, dt},
                                  {0, 0, 1}});
constexpr numeric_matrix<3, 1> B;
constexpr numeric_matrix<1, 3> H({{1}});
constexpr numeric_matrix<3, 3> Q = numeric_matrix<3, 3>::diagonals(base_noise_value);
constexpr numeric_matrix<1, 1> R = numeric_matrix<1, 1>::diagonals(base_noise_value);
constexpr numeric_vector<3> x0;// {x, v, a}

kalman_filter_t<3, 1, 1> kf(F, B, H, Q, R, x0);

int main() {
    constexpr size_t SAMPLES                = 16;
    numeric_vector<1> controls[SAMPLES]     = {};
    numeric_vector<1> measurements[SAMPLES] = {};

    real_t i_r = 0;
    for (auto &x: measurements) x = {i_r++};

    for (size_t i = 0; i < SAMPLES; ++i) {
        kf.predict(controls[i]);
        std::cout << "Predicted state: ";
        for (auto &state: kf.state_vector) {
            std::cout << state << " ";
        }
        std::cout << '\n';

        kf.update(measurements[i]);
        std::cout << i + 1 << " - Estimated state: ";
        for (auto &state: kf.state_vector) {
            std::cout << state << " ";
        }
        std::cout << '\n';
    }

    return 0;
}
