#include "vt_linalg"
#include "vt_kalman"

using namespace vt;

real_t dt = 1.0;
real_t hdts = 0.5 * dt * dt;
real_t base_noise_value = 0.001;

numeric_matrix<9, 9> F({{1, 0, 0, dt, 0,  0,  hdts, 0,    0},
                        {0, 1, 0, 0,  dt, 0,  0,    hdts, 0},
                        {0, 0, 1, 0,  0,  dt, 0,    0,    hdts},
                        {0, 0, 0, 1,  0,  0,  dt,   0,    0},
                        {0, 0, 0, 0,  1,  0,  0,    dt,   0},
                        {0, 0, 0, 0,  0,  1,  0,    0,    dt},
                        {0, 0, 0, 0,  0,  0,  1,    0,    0},
                        {0, 0, 0, 0,  0,  0,  0,    1,    0},
                        {0, 0, 0, 0,  0,  0,  0,    0,    1}});
numeric_matrix<9, 1> B;
numeric_matrix<3, 9> H({{1},
                        {0, 1},
                        {0, 0, 1}});
numeric_matrix<9, 9> Q = numeric_matrix<9, 9>::diagonals(base_noise_value);
numeric_matrix<3, 3> R = numeric_matrix<3, 3>::diagonals(base_noise_value);
numeric_vector<9> x0; // {x, y, z, v_x, v_y, v_z, a_x, a_y, a_z}

kalman_filter_t<9, 3, 1> kf(F, B, H, Q, R, x0);

int main() {
    return 0;
}
