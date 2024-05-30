#ifndef VT_LINALG_KALMAN_LUT_H
#define VT_LINALG_KALMAN_LUT_H

#include "kalman.h"

namespace vt {
    class kf_basic_1d : public kalman_filter_t<3, 1, 1> {
    public:
        kf_basic_1d(real_t dt, real_t base_noise_value)
            : kalman_filter_t<3, 1, 1>(
                      make_numeric_matrix<3, 3>({{1, dt, 0.5 * dt * dt},
                                                 {0, 1, dt},
                                                 {0, 0, 1}}),
                      make_numeric_matrix<3, 1>(),
                      make_numeric_matrix<1, 3>({{1}}),
                      numeric_matrix<3, 3>::diagonals(base_noise_value),
                      numeric_matrix<1, 1>::diagonals(base_noise_value),
                      make_numeric_vector<3>()) {}
    };

    class kf_basic_2d : public kalman_filter_t<6, 2, 1> {
    public:
        kf_basic_2d(real_t dt, real_t base_noise_value)
            : kalman_filter_t<6, 2, 1>(
                      make_numeric_matrix<6, 6>({{1, 0, dt, 0, 0.5 * dt * dt, 0},
                                                 {0, 1, 0, dt, 0, 0.5 * dt * dt},
                                                 {0, 0, 1, 0, dt, 0},
                                                 {0, 0, 0, 1, 0, dt},
                                                 {0, 0, 0, 0, 1, 0},
                                                 {0, 0, 0, 0, 0, 1}}),
                      make_numeric_matrix<6, 1>(),
                      make_numeric_matrix<2, 6>({{1},
                                                 {0, 1}}),
                      numeric_matrix<6, 6>::diagonals(base_noise_value),
                      numeric_matrix<2, 2>::diagonals(base_noise_value),
                      make_numeric_vector<6>()) {}
    };

    class kf_basic_3d : public kalman_filter_t<9, 3, 1> {
    public:
        kf_basic_3d(real_t dt, real_t base_noise_value)
            : kalman_filter_t<9, 3, 1>(
                      make_numeric_matrix<9, 9>({{1, 0, 0, dt, 0, 0, 0.5 * dt * dt, 0, 0},
                                                 {0, 1, 0, 0, dt, 0, 0, 0.5 * dt * dt, 0},
                                                 {0, 0, 1, 0, 0, dt, 0, 0, 0.5 * dt * dt},
                                                 {0, 0, 0, 1, 0, 0, dt, 0, 0},
                                                 {0, 0, 0, 0, 1, 0, 0, dt, 0},
                                                 {0, 0, 0, 0, 0, 1, 0, 0, dt},
                                                 {0, 0, 0, 0, 0, 0, 1, 0, 0},
                                                 {0, 0, 0, 0, 0, 0, 0, 1, 0},
                                                 {0, 0, 0, 0, 0, 0, 0, 0, 1}}),
                      make_numeric_matrix<9, 1>(),
                      make_numeric_matrix<3, 9>({{1},
                                                 {0, 1},
                                                 {0, 0, 1}}),
                      numeric_matrix<9, 9>::diagonals(base_noise_value),
                      numeric_matrix<3, 3>::diagonals(base_noise_value),
                      make_numeric_vector<9>()) {}
    };
}  // namespace vt

#endif  //VT_LINALG_KALMAN_LUT_H
