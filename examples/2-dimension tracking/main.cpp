#include "vt_linalg"
#include "vt_kalman"

using namespace vt;

real_t dt = 1.0;
real_t hdts = 0.5 * dt * dt;
real_t base_noise_value = 0.001;

numeric_matrix<6, 6> F({{1, 0, dt, 0,  hdts, 0},
                        {0, 1, 0,  dt, 0,    hdts},
                        {0, 0, 1,  0,  dt,   0},
                        {0, 0, 0,  1,  0,    dt},
                        {0, 0, 0,  0,  1,    0},
                        {0, 0, 0,  0,  0,    1}});
numeric_matrix<6, 1> B;
numeric_matrix<2, 6> H({{0, 1}});
numeric_matrix<6, 6> Q = numeric_matrix<6, 6>::diagonals(base_noise_value);
numeric_matrix<2, 2> R = numeric_matrix<2, 2>::diagonals(base_noise_value);
numeric_vector<6> x0; // {x, y, v_x, v_y, a_x, a_y}

kalman_filter_t<6, 2, 1> kf(F, B, H, Q, R, x0);

int main() {
    return 0;
}
