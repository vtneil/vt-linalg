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
            size_t MeasurementVectorDimension,
            size_t ControlVectorDimension>
    class simple_kalman_filter_t {
    private:
        // Note that numeric_matrix<N_, M_> maps from R_^M_ to R_^N_
        static constexpr size_t N_ = StateVectorDimension;       // ALias
        static constexpr size_t M_ = MeasurementVectorDimension; // Alias
        static constexpr size_t L_ = ControlVectorDimension;     // Alias

    protected:
        numeric_matrix<N_, N_> &F_; // state-transition model
        numeric_matrix<N_, L_> &B_; // control-input model
        numeric_matrix<M_, N_> &H_; // measurement model
        numeric_matrix<N_, N_> &Q_; // covariance of the process noise
        numeric_matrix<M_, M_> &R_; // covariance of the measurement noise
        numeric_vector<N_> &x_;    // state vector
        numeric_matrix<N_, N_> P_; // state covariance, self-initialized as Q_

    public:
        /**
         * Simple Kalman filter array-copying constructor
         *
         * @param F_matrix state-transition model
         * @param B_matrix control-input model
         * @param H_matrix measurement model
         * @param Q_matrix covariance of the process noise
         * @param R_matrix covariance of the measurement noise
         * @param x_0 initial state vector
         */
        simple_kalman_filter_t(
                numeric_matrix<N_, N_> &F_matrix,
                numeric_matrix<N_, L_> &B_matrix,
                numeric_matrix<M_, N_> &H_matrix,
                numeric_matrix<N_, N_> &Q_matrix,
                numeric_matrix<M_, M_> &R_matrix,
                numeric_vector<N_> &x_0
        ) : F_(F_matrix), B_(B_matrix), H_(H_matrix), Q_(Q_matrix), R_(R_matrix), x_(x_0), P_(Q_matrix) {}

        simple_kalman_filter_t(const simple_kalman_filter_t &other) = default;

        simple_kalman_filter_t(simple_kalman_filter_t &&other) noexcept = default;

        /**
         * Kalman filter prediction
         *
         * @param u control input vector
         */
        void predict(const numeric_vector<L_> &u) {
            x_ = F_ * x_ + B_ * u;
            P_ = F_ * P_ * F_.transpose() + Q_;
        }

        /**
         * Kalman filter update
         *
         * @param z Measurement vector
         */
        void update(const numeric_vector<M_> &z) {
            numeric_vector<M_> y = z - H_ * x_;
            numeric_matrix<N_, M_> K = P_ * H_.transpose() * (H_ * P_ * H_.transpose() + R_).inverse();
            x_ = x_ + K * y;
            P_ = (numeric_matrix<N_, N_>::identity() - K * H_) * P_;
        }

        __VT_FORCE_INLINE numeric_vector<N_> &state_vector() { return x_; }
    };

    template<size_t StateVectorDimension,
            size_t MeasurementVectorDimension,
            size_t ControlVectorDimension>
    using kalman_filter_t =
            simple_kalman_filter_t<StateVectorDimension, MeasurementVectorDimension, ControlVectorDimension>;
}

#endif //VNET_LINALG_KALMAN_H
