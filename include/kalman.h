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

#ifndef VT_LINALG_KALMAN_H
#define VT_LINALG_KALMAN_H

#include "numeric_matrix.h"
#include "numeric_vector.h"
#include "standard_utility.h"

namespace vt {
    template<size_t StateVectorDimension, size_t MeasurementVectorDimension, size_t ControlVectorDimension>
    class simple_kalman_filter_t {
    private:
        // Note that numeric_matrix<N_, M_> maps from R_^M_ to R_^N_
        static constexpr size_t N_ = StateVectorDimension;        // ALias
        static constexpr size_t M_ = MeasurementVectorDimension;  // Alias
        static constexpr size_t L_ = ControlVectorDimension;      // Alias

    protected:
        const numeric_matrix<N_, N_> &F_;  // state-transition model
        const numeric_matrix<N_, L_> &B_;  // control-input model
        const numeric_matrix<M_, N_> &H_;  // measurement model
        const numeric_matrix<N_, N_> &Q_;  // covariance of the process noise
        const numeric_matrix<M_, M_> &R_;  // covariance of the measurement noise
        numeric_vector<N_> x_;             // state vector
        numeric_matrix<N_, N_> P_;         // state covariance, self-initialized as Q_

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
        constexpr simple_kalman_filter_t(
                const numeric_matrix<N_, N_> &F_matrix,
                const numeric_matrix<N_, L_> &B_matrix,
                const numeric_matrix<M_, N_> &H_matrix,
                const numeric_matrix<N_, N_> &Q_matrix,
                const numeric_matrix<M_, M_> &R_matrix,
                const numeric_vector<N_> &x_0)
            : F_{F_matrix}, B_{B_matrix}, H_{H_matrix},
              Q_{Q_matrix}, R_{R_matrix}, x_{x_0}, P_{Q_matrix} {}

        constexpr simple_kalman_filter_t(const simple_kalman_filter_t &) = delete;

        constexpr simple_kalman_filter_t(simple_kalman_filter_t &&) noexcept = delete;

        /**
         * Kalman filter prediction
         *
         * @param u control input vector
         */
        simple_kalman_filter_t &predict(const numeric_vector<L_> &u = {}) {
            x_ = vt::move(F_ * x_ + B_ * u);
            P_ = vt::move(F_ * P_.matmul_T(F_) + Q_);
            return *this;
        }

        /**
         * Kalman filter update
         *
         * @param z Measurement vector
         */
        simple_kalman_filter_t &update(const numeric_vector<M_> &z) {
            numeric_vector<M_> y_        = vt::move(z - H_ * x_);
            numeric_matrix<N_, M_> P_H_t = vt::move(P_.matmul_T(H_));
            numeric_matrix<M_, M_> S_    = vt::move(H_ * P_H_t + R_);
            numeric_matrix<N_, M_> K_    = vt::move(P_H_t * S_.inverse());

            x_ += K_ * y_;
            P_ = vt::move((numeric_matrix<N_, N_>::identity() - K_ * H_) * P_);

            return *this;
        }

        template<typename... Ts>
        simple_kalman_filter_t &update(Ts... vs) { return update(make_numeric_vector({vs...})); }

        const numeric_vector<N_> &state_vector = x_;

        const real_t &state = x_[0];
    };

    template<size_t StateVectorDimension, size_t MeasurementVectorDimension, size_t ControlVectorDimension>
    class adaptive_kalman_filter_t {
    private:
        // Note that numeric_matrix<N_, M_> maps from R_^M_ to R_^N_
        static constexpr size_t N_ = StateVectorDimension;        // ALias
        static constexpr size_t M_ = MeasurementVectorDimension;  // Alias
        static constexpr size_t L_ = ControlVectorDimension;      // Alias

    protected:
        const numeric_matrix<N_, N_> &F_;  // state-transition model
        const numeric_matrix<N_, L_> &B_;  // control-input model
        const numeric_matrix<M_, N_> &H_;  // measurement model
        numeric_matrix<N_, N_> &Q_;        // covariance of the process noise
        numeric_matrix<M_, M_> &R_;        // covariance of the measurement noise
        numeric_vector<N_> x_;             // state vector
        numeric_matrix<N_, N_> P_;         // state covariance, self-initialized as Q_

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
        constexpr adaptive_kalman_filter_t(
                const numeric_matrix<N_, N_> &F_matrix,
                const numeric_matrix<N_, L_> &B_matrix,
                const numeric_matrix<M_, N_> &H_matrix,
                numeric_matrix<N_, N_> &Q_matrix,
                numeric_matrix<M_, M_> &R_matrix,
                const numeric_vector<N_> &x_0)
            : F_{F_matrix}, B_{B_matrix}, H_{H_matrix},
              Q_{Q_matrix}, R_{R_matrix}, x_{x_0}, P_{Q_matrix} {}

        constexpr adaptive_kalman_filter_t(const adaptive_kalman_filter_t &) = delete;

        constexpr adaptive_kalman_filter_t(adaptive_kalman_filter_t &&) noexcept = delete;

        /**
         * Kalman filter prediction
         *
         * @param u control input vector
         */
        adaptive_kalman_filter_t &predict(const numeric_vector<L_> &u = {}) {
            x_ = vt::move(F_ * x_ + B_ * u);
            P_ = vt::move(F_ * P_.matmul_T(F_) + Q_);
            return *this;
        }

        /**
         * Kalman filter update
         *
         * @param z Measurement vector
         */
        adaptive_kalman_filter_t &update(const numeric_vector<M_> &z) {
            numeric_vector<M_> y_        = vt::move(z - H_ * x_);
            numeric_matrix<N_, M_> P_H_t = vt::move(P_.matmul_T(H_));
            numeric_matrix<M_, M_> S_    = vt::move(H_ * P_H_t + R_);
            numeric_matrix<N_, M_> K_    = vt::move(P_H_t * S_.inverse());

            x_ += K_ * y_;
            P_ = vt::move((numeric_matrix<N_, N_>::identity() - K_ * H_) * P_);

            numeric_matrix<M_, 1> y_mat = y_.as_matrix_col();
            numeric_matrix<M_, M_> y_yT = y_mat.matmul_T(y_mat);

            adapt_R(y_yT, S_);
            adapt_Q(K_, y_yT);

            return *this;
        }

        template<typename... Ts>
        adaptive_kalman_filter_t &update(Ts... vs) { return update(make_numeric_vector({vs...})); }

        const numeric_vector<N_> &state_vector = x_;

        const real_t &state = x_[0];

    private:
        void adapt_R(const numeric_matrix<M_, M_> &y_yT, const numeric_matrix<M_, M_> &S) {
            static constexpr real_t alpha = 0.1;  // Smoothing factor
            R_                            = (1 - alpha) * R_ + alpha * (y_yT + S);
        }

        void adapt_Q(const numeric_matrix<N_, M_> &K, const numeric_matrix<M_, M_> &y_yT) {
            static constexpr real_t beta = 0.1;  // Smoothing factor
            Q_                           = (1 - beta) * Q_ + beta * (K * y_yT * K.transpose());
        }
    };

    template<size_t StateVectorDimension, size_t MeasurementVectorDimension, size_t ControlVectorDimension>
    class extended_kalman_filter_t {
    private:
        // Note that numeric_matrix<N_, M_> maps from R_^M_ to R_^N_
        static constexpr size_t N_ = StateVectorDimension;        // ALias
        static constexpr size_t M_ = MeasurementVectorDimension;  // Alias
        static constexpr size_t L_ = ControlVectorDimension;      // Alias

    public:
        typedef numeric_vector<N_> (*state_func_t)(const numeric_vector<N_> &x, const numeric_vector<L_> &u);

        typedef numeric_matrix<N_, N_> (*state_jacobian_t)(const numeric_vector<N_> &x, const numeric_vector<L_> &u);

        typedef numeric_vector<M_> (*observation_func_t)(const numeric_vector<N_> &x);

        typedef numeric_matrix<M_, N_> (*observation_jacobian_t)(const numeric_vector<N_> &x);

    protected:
        const state_func_t f_;             // state-transition model
        const state_jacobian_t Fj_;        // state-transition Jacobian
        const observation_func_t h_;       // measurement model
        const observation_jacobian_t Hj_;  // measurement Jacobian
        const numeric_matrix<N_, N_> &Q_;  // covariance of the process noise
        const numeric_matrix<M_, M_> &R_;  // covariance of the measurement noise
        numeric_vector<N_> x_;             // state vector
        numeric_matrix<N_, N_> P_;         // state covariance, self-initialized as Q_

    public:
        constexpr extended_kalman_filter_t(
                const state_func_t f_vec_func,
                const state_jacobian_t Fj_mat_func,
                const observation_func_t h_vec_func,
                const observation_jacobian_t Hj_mat_func,
                const numeric_matrix<N_, N_> &Q_matrix,
                const numeric_matrix<M_, M_> &R_matrix,
                const numeric_vector<N_> &x_0)
            : f_(f_vec_func), Fj_{Fj_mat_func}, h_{h_vec_func}, Hj_{Hj_mat_func},
              Q_{Q_matrix}, R_{R_matrix}, x_{x_0}, P_{Q_matrix} {}

        /**
         * Disabled copy constructor
         */
        constexpr extended_kalman_filter_t(const extended_kalman_filter_t &) = delete;

        /**
         * Disabled move constructor
         */
        constexpr extended_kalman_filter_t(extended_kalman_filter_t &&) noexcept = delete;

        extended_kalman_filter_t &predict(const numeric_vector<L_> &u = {}) {
            x_                        = vt::move(f_(x_, u));
            numeric_matrix<N_, N_> F_ = vt::move(Fj_(x_, u));
            P_                        = vt::move(F_ * P_.matmul_T(F_) + Q_);
            return *this;
        }

        extended_kalman_filter_t &update(const numeric_vector<M_> &z) {
            numeric_vector<M_> y_          = vt::move(z - h_(x_));
            numeric_matrix<M_, N_> Hjx_    = vt::move(Hj_(x_));
            numeric_matrix<N_, M_> P_Hjx_t = vt::move(P_.matmul_T(Hjx_));
            numeric_matrix<M_, M_> S_      = vt::move(Hjx_ * P_Hjx_t + R_);
            numeric_matrix<N_, M_> K_      = vt::move(P_Hjx_t * S_.inverse());
            x_ += K_ * y_;
            P_ = vt::move((numeric_matrix<N_, N_>::identity() - K_ * Hj_(x_)) * P_);
            return *this;
        }

        template<typename... Ts>
        extended_kalman_filter_t &update(Ts... vs) { return update(make_numeric_vector({vs...})); }

        const numeric_vector<N_> &state_vector = x_;

        const real_t &state = x_[0];
    };

    namespace future {
        template<size_t StateVectorDimension, size_t MeasurementVectorDimension, size_t ControlVectorDimension>
        class unscented_kalman_filter_t {
        private:
            // Note that numeric_matrix<N_, M_> maps from R_^M_ to R_^N_
            static constexpr size_t N_ = StateVectorDimension;          // ALias
            static constexpr size_t M_ = MeasurementVectorDimension;    // Alias
            static constexpr size_t L_ = ControlVectorDimension;        // Alias
            static constexpr size_t Z_ = 2 * StateVectorDimension + 1;  // Sigma Points
        public:
            typedef numeric_vector<N_> (*state_func_t)(const numeric_vector<N_> &x, const numeric_vector<L_> &u);

            typedef numeric_vector<M_> (*observation_func_t)(const numeric_vector<N_> &x);

        protected:
            const state_func_t f_;             // state-transition model
            const observation_func_t h_;       // measurement model
            const numeric_matrix<N_, N_> &Q_;  // covariance of the process noise
            const numeric_matrix<M_, M_> &R_;  // covariance of the measurement noise
            numeric_vector<N_> x_;             // state vector
            numeric_matrix<N_, N_> P_;         // state covariance, self-initialized as Q_

        public:
            void predict() {
            }

            void update() {
            }

        protected:
            void cp_sigma() {
            }

            void cp_weight() {
            }
        };
    }  // namespace future

    // Aliases

    template<size_t StateVectorDimension, size_t MeasurementVectorDimension, size_t ControlVectorDimension>
    using kalman_filter_t =
            simple_kalman_filter_t<StateVectorDimension, MeasurementVectorDimension, ControlVectorDimension>;

    template<size_t StateVectorDimension, size_t MeasurementVectorDimension, size_t ControlVectorDimension>
    using KF =
            simple_kalman_filter_t<StateVectorDimension, MeasurementVectorDimension, ControlVectorDimension>;

    template<size_t StateVectorDimension, size_t MeasurementVectorDimension, size_t ControlVectorDimension>
    using EKF =
            extended_kalman_filter_t<StateVectorDimension, MeasurementVectorDimension, ControlVectorDimension>;
}  // namespace vt

#endif  //VT_LINALG_KALMAN_H
