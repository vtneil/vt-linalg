/**
 * @file kalman.h
 * @author Vivatsathorn Thitasirivit
 * @date 19 June 2023
 * @brief Kalman filter implementation library
 *
 * For Simple KF, using Discrete-time KF.\n
 *
 * Referencing https://www.intechopen.com/chapters/63164
 */

#ifndef VNET_LINALG_KALMAN_H
#define VNET_LINALG_KALMAN_H

#include "standard_utility.h"
#include "numeric_vector.h"
#include "numeric_matrix.h"

namespace vt {
    template<size_t StateVectorDimension,
            size_t ObservationVectorDimension,
            size_t ControlVectorDimension>
    class simple_kalman_filter_t {
    private:
        // Note that numeric_matrix<N, M> maps from R^M to R^N
        static constexpr size_t N = StateVectorDimension;       // ALias
        static constexpr size_t M = ObservationVectorDimension; // Alias
        static constexpr size_t L = ControlVectorDimension;     // Alias

    protected:
        numeric_matrix<N, N> &F; // state-transition model
        numeric_matrix<N, L> &B; // control-input model
        numeric_matrix<M, N> &H; // measurement model
        numeric_matrix<N, N> &Q; // covariance of the process noise
        numeric_matrix<M, M> &R; // covariance of the measurement noise
        numeric_vector<N> &x;    // state vector
        numeric_matrix<N, N> P; // state covariance, self-initialized as Q

    public:
        /**
         * Simple Kalman filter array-copying constructor
         *
         * @param F state-transition model
         * @param B control-input model
         * @param H measurement model
         * @param Q covariance of the process noise
         * @param R covariance of the measurement noise
         * @param x_0 initial state vector
         */
        simple_kalman_filter_t(
                numeric_matrix<N, N> &F,
                numeric_matrix<N, L> &B,
                numeric_matrix<M, N> &H,
                numeric_matrix<N, N> &Q,
                numeric_matrix<M, M> &R,
                numeric_vector<N> &x_0
        ) : F(F), B(B), H(H), Q(Q), R(R), x(x_0), P(Q) {}

        simple_kalman_filter_t(const simple_kalman_filter_t &other) = default;

        simple_kalman_filter_t(simple_kalman_filter_t &&other) noexcept = default;

        /**
         * Kalman filter prediction
         *
         * @param u control input vector
         */
        void predict(const numeric_vector<L> &u) {
            x = F * x + B * u;
            P = F * P * F.transpose() + Q;
        }

        /**
         * Kalman filter update
         *
         * @param z Measurement vector
         */
        void update(const numeric_vector<M> &z) {
            numeric_vector<M> y = z - H * x;
            numeric_matrix<N, M> K = P * H.transpose() * (H * P * H.transpose() + R).inverse();
            x = x + K * y;
            P = (numeric_matrix<N, N>::identity() - K * H) * P;
        }

        __VT_FORCE_INLINE numeric_vector<N> &state_vector() { return x; }
    };

    template<size_t StateVectorDimension,
            size_t ObservationVectorDimension,
            size_t ControlVectorDimension>
    using kalman_filter_t =
            simple_kalman_filter_t<StateVectorDimension, ObservationVectorDimension, ControlVectorDimension>;
}

#endif //VNET_LINALG_KALMAN_H
