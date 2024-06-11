#ifndef VT_LINALG_KALMAN_LUT_H
#define VT_LINALG_KALMAN_LUT_H

#include "kalman.h"

namespace vt {
    /**
     * Variable dt wrapper for Kalman filter
     *
     * @tparam Degree Max degree to calculate d^N/dt^N
     */
    template<size_t Degree>
    class vdt {
    protected:
        real_t m_dt[Degree] = {};

    public:
        vdt() = delete;

        vdt(const vdt &other) = default;

        vdt(vdt &&other) noexcept = default;

        explicit vdt(const real_t &dt) { update_dt(dt); }

        void update_dt(const real_t &new_dt) { update_dt_helper<0>(new_dt); }

        numeric_matrix<Degree + 1, Degree + 1> generate_F() {
            numeric_matrix<Degree + 1, Degree + 1> F_out = {};

            for (size_t i = 0; i < Degree + 1; ++i) {
                for (size_t j = 0; j < Degree + 1; ++j) {
                    if (i == j) {
                        F_out[i][j] = 1.;
                    } else if (i < j) {
                        F_out[i][j] = j - i <= Degree
                                              ? m_dt[j - i - 1]
                                              : 0.;
                    } else {
                        F_out[i][j] = 0.;
                    }
                }
            }

            return F_out;
        }

    protected:
        template<size_t Index>
        void update_dt_helper(const real_t &new_dt) {
            if constexpr (Index == 0) {
                m_dt[0] = new_dt;
                update_dt_helper<1>(new_dt);
            } else if constexpr (Index == 1) {
                m_dt[1] = 0.5 * m_dt[0] * m_dt[0];
                update_dt_helper<2>(new_dt);
            } else if constexpr (Index < Degree) {
                m_dt[Index] = m_dt[Index - 1] * m_dt[Index - 2] / static_cast<real_t>(Index + 1);
                update_dt_helper<Index + 1>(new_dt);
            }
        }
    };

    template<size_t Order>
    class kf_pos {
    public:
        static_assert(Order > 0, "Order must be non-zero.");

        using vdt_type = vdt<Order - 1>;
        using kf_type  = adaptive_kalman_filter_t<Order, 1, 1>;

        vdt_type vdt;
        numeric_matrix<Order, Order> F;
        numeric_matrix<Order, 1> B;
        numeric_matrix<1, Order> H;
        numeric_matrix<Order, Order> Q;
        numeric_matrix<1, 1> R;
        kf_type kf;

        kf_pos(const real_t &dt, const real_t &covariance,
               const real_t &alpha, const real_t &beta)
            : vdt{vdt_type(dt)},
              F{vdt.generate_F()},
              B{make_numeric_matrix<Order, 1>()},
              H{make_numeric_matrix<1, Order>({{1}})},
              Q{numeric_matrix<Order, Order>::diagonals(covariance)},
              R{numeric_matrix<1, 1>::diagonals(covariance)},
              kf{kf_type(F, B, H, Q, R, make_numeric_vector<Order>(), alpha, beta)} {}

        kf_pos(const kf_pos &) = default;

        kf_pos(kf_pos &&) noexcept = default;

        kf_pos &operator=(const kf_pos &) = default;

        kf_pos &operator=(kf_pos &&) noexcept = default;

        void update_dt(const real_t &new_dt) {
            vdt.update_dt(new_dt);
            F = vdt.generate_F();
        }
    };

    template<size_t Order>
    class kf_acc {
    public:
        static_assert(Order >= 3, "Order must be at least 3 to include acceleration.");

        using vdt_type = vdt<Order - 1>;
        using kf_type  = adaptive_kalman_filter_t<Order, 1, 1>;

        vdt_type vdt;
        numeric_matrix<Order, Order> F;
        numeric_matrix<Order, 1> B;
        numeric_matrix<1, Order> H;
        numeric_matrix<Order, Order> Q;
        numeric_matrix<1, 1> R;
        kf_type kf;

        kf_acc(const real_t &dt, const real_t &covariance,
               const real_t &alpha, const real_t &beta)
            : vdt{vdt_type(dt)},
              F{vdt.generate_F()},
              B{make_numeric_matrix<Order, 1>()},
              H{make_numeric_matrix<1, Order>({{0, 0, 1}})},
              Q{numeric_matrix<Order, Order>::diagonals(covariance)},
              R{numeric_matrix<1, 1>::diagonals(covariance)},
              kf{kf_type(F, B, H, Q, R, make_numeric_vector<Order>(), alpha, beta)} {}

        kf_acc(const kf_acc &) = default;

        kf_acc(kf_acc &&) noexcept = default;

        kf_acc &operator=(const kf_acc &) = default;

        kf_acc &operator=(kf_acc &&) noexcept = default;

        void update_dt(const real_t &new_dt) {
            vdt.update_dt(new_dt);
            F = vdt.generate_F();
        }
    };

    template<size_t Order>
    class kf_pos_acc {
    public:
        static_assert(Order >= 3, "Order must be at least 3 to include acceleration.");

        using vdt_type = vdt<3>;
        using kf_type  = adaptive_kalman_filter_t<4, 2, 1>;

        vdt_type vdt;
        numeric_matrix<Order, Order> F;
        numeric_matrix<Order, 1> B;
        numeric_matrix<2, Order> H;
        numeric_matrix<Order, Order> Q;
        numeric_matrix<2, 2> R;
        kf_type kf;

        kf_pos_acc(const real_t &dt, const real_t &covariance,
                   const real_t &alpha, const real_t &beta)
            : vdt{vdt_type(dt)},
              F{vdt.generate_F()},
              B{make_numeric_matrix<Order, 1>()},
              H{make_numeric_matrix<2, Order>({{1, 0, 0}, {0, 0, 1}})},
              Q{numeric_matrix<Order, Order>::diagonals(covariance)},
              R{numeric_matrix<2, 2>::diagonals(covariance)},
              kf{kf_type(F, B, H, Q, R, make_numeric_vector<Order>(), alpha, beta)} {}

        kf_pos_acc(const kf_pos_acc &) = default;

        kf_pos_acc(kf_pos_acc &&) noexcept = default;

        kf_pos_acc &operator=(const kf_pos_acc &) = default;

        kf_pos_acc &operator=(kf_pos_acc &&) noexcept = default;

        void update_dt(const real_t &new_dt) {
            vdt.update_dt(new_dt);
            F = vdt.generate_F();
        }
    };
}  // namespace vt

#endif  //VT_LINALG_KALMAN_LUT_H
