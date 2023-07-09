/**
 * Jerk (3rd integral) simple kalman filter model.
 */

#include "vt_linalg"
#include "vt_kalman"

using namespace vt;

constexpr real_t dt = 1.0;
constexpr real_t d2t = integral_coefficient<1>() * (dt * dt);
constexpr real_t d3t = integral_coefficient<2>() * (dt * dt * dt);
constexpr real_t base_noise_value = 0.001;

constexpr numeric_matrix<4, 4> F = make_numeric_matrix({{1, dt, d2t, d3t},
                                                        {0, 1,  dt,  d2t},
                                                        {0, 0,  1,   dt},
                                                        {0, 0,  0,   1}});
constexpr numeric_matrix<4, 1> B;
constexpr numeric_matrix<1, 4> H = make_numeric_matrix({{1, 0, 0, 0}});
constexpr numeric_matrix<4, 4> Q = numeric_matrix<4, 4>::diagonals(base_noise_value);
constexpr numeric_matrix<1, 1> R = numeric_matrix<1, 1>::diagonals(base_noise_value);
constexpr numeric_vector<4> x0; // {x, v, a, j}

kalman_filter_t<4, 1, 1> kf(F, B, H, Q, R, x0);

int main() {
    return 0;
}
