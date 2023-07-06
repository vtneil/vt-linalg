#ifndef VT_LINALG_KALMAN_LUT_H
#define VT_LINALG_KALMAN_LUT_H

namespace vt {
    namespace kalman {
        constexpr numeric_matrix<3, 3> make_F_1D(const real_t dt) {
            return numeric_matrix < 3, 3 > (
                    {
                        { 1, dt, (0.5 * dt * dt) },
                        { 0, 1, dt },
                        { 0, 0, 1 }
                    }
            );
        }
    }
}

#endif //VT_LINALG_KALMAN_LUT_H
