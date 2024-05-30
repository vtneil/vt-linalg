#include "vt_linalg"
#include "vt_kalman"

using namespace vt;

real_t dt1 = 1.0;
real_t dt2 = 0.5 * dt1 * dt1;
real_t base_noise_value = 0.001;

numeric_matrix<3, 3> F({{1, dt1, dt2},
                        {0, 1,  dt1},
                        {0, 0,  1}});
numeric_matrix<3, 1> B;
numeric_matrix<1, 3> H({{1}});
numeric_matrix<3, 3> Q = numeric_matrix<3, 3>::diagonals(base_noise_value);
numeric_matrix<1, 1> R = numeric_matrix<1, 1>::diagonals(base_noise_value);
numeric_vector<3> x0; // {x, v, a}

kalman_filter_t<3, 1, 1> kf(F, B, H, Q, R, x0);

int main() {
    return 0;
}
