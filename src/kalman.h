/**
 * @file kalman.h
 * @author Vivatsathorn Thitasirivit
 * @date 19 June 2023
 * @brief Kalman filter implementation library
 *
 * For Simple KF, using Discrete-time KF.\n
 *
 * The equations are:\n
 * x[k+1] = A.x[k] + B.u[k] + G.w[k]\n
 * y[k] = C.x[k] + D.u[k] + v[k]\n
 *
 * w is noise vector ~ N(0, Q).\n
 * v is noise vector ~ N(0, R).\n
 * x is state vector.\n
 * y is measurement vector.\n
 * A is state-transition matrix.\n
 * B is control matrix for state.\n
 * G is noise matrix.\n
 * C is measurement matrix.\n
 * D is control matrix for measurement.\n
 */

#ifndef VNET_LINALG_KALMAN_H
#define VNET_LINALG_KALMAN_H

#include "standard_utility.h"
#include "numeric_vector.h"
#include "numeric_matrix.h"

namespace vt {
    template<size_t StateVectorDimension, size_t ObservationVectorDimension>
    class simple_kalman_filter_t {
    private:
        // Note that numeric_matrix<N, M> maps from R^M to R^N
        static constexpr size_t N = StateVectorDimension;       // ALias
        static constexpr size_t M = ObservationVectorDimension; // Alias

    protected:
        numeric_matrix<N, N> F; // state-transition model
        numeric_matrix<M, N> H; // observation model
        numeric_matrix<N, N> Q; // covariance of the process noise
        numeric_matrix<M, M> R; // covariance of the observation noise
        numeric_vector<N> x;    // estimated state
        numeric_matrix<N, N> P; // state covariance

    public:
        simple_kalman_filter_t(
                const numeric_matrix<N, N> &F,
                const numeric_matrix<M, N> &H,
                const numeric_matrix<N, N> &Q,
                const numeric_matrix<M, M> &R,
                const numeric_vector<N> &x0
        ) : F(F), H(H), Q(Q), R(R), x(x0), P(Q) {}

        void predict() {
            x = F * x;
            P = F * P * F.transpose() + Q;
        }

        void update(const numeric_vector<M> &z) {
            numeric_vector<M> y = z - H * x;
            numeric_matrix<M, M> S = H * P * H.transpose() + R;
            numeric_matrix<N, M> K = P * H.transpose() * S.inv();
            x = x + K * y;
            P = (numeric_matrix<N, N>::identity() - K * H) * P;
        }

        numeric_vector<N> &state_vector() { return x; }
    };
}

#endif //VNET_LINALG_KALMAN_H
