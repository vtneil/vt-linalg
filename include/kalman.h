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
    class kalman_filter_t {
    private:
        // Note that numeric_matrix<N_, M_> maps from R_^M_ to R_^N_
        static constexpr size_t N_ = StateVectorDimension;        // ALias
        static constexpr size_t M_ = MeasurementVectorDimension;  // Alias
        static constexpr size_t L_ = ControlVectorDimension;      // Alias

    protected:
        numeric_matrix<N_, N_> &F_;        // state-transition model
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
        constexpr kalman_filter_t(
                numeric_matrix<N_, N_> &F_matrix,
                const numeric_matrix<N_, L_> &B_matrix,
                const numeric_matrix<M_, N_> &H_matrix,
                const numeric_matrix<N_, N_> &Q_matrix,
                const numeric_matrix<M_, M_> &R_matrix,
                const numeric_vector<N_> &x_0,
                const real_t & = 0.,
                const real_t & = 0.)
            : F_{F_matrix}, B_{B_matrix}, H_{H_matrix},
              Q_{Q_matrix}, R_{R_matrix}, x_{x_0}, P_{Q_matrix} {}

        constexpr kalman_filter_t(const kalman_filter_t &) = default;

        constexpr kalman_filter_t(kalman_filter_t &&) noexcept = default;

        kalman_filter_t &operator=(const kalman_filter_t &) = default;

        kalman_filter_t &operator=(kalman_filter_t &&) noexcept = default;

        /**
         * Kalman filter prediction
         *
         * @param u control input vector
         */
        kalman_filter_t &predict(const numeric_vector<L_> &u = {}) {
            x_ = vt::move(F_ * x_ + B_ * u);
            P_ = vt::move(F_ * P_.matmul_T(F_) + Q_);
            return *this;
        }

        /**
         * Kalman filter update
         *
         * @param z Measurement vector
         */
        kalman_filter_t &update(const numeric_vector<M_> &z) {
            numeric_vector<M_> y_        = vt::move(z - H_ * x_);
            numeric_matrix<N_, M_> P_H_t = vt::move(P_.matmul_T(H_));
            numeric_matrix<M_, M_> S_    = vt::move(H_ * P_H_t + R_);
            numeric_matrix<N_, M_> K_    = vt::move(P_H_t * S_.inverse());

            x_ += K_ * y_;
            P_ = vt::move((numeric_matrix<N_, N_>::identity() - K_ * H_) * P_);

            return *this;
        }

        kalman_filter_t &operator<<(const numeric_vector<M_> &z) {
            return predict().update(z);
        }

        template<typename... Ts>
        kalman_filter_t &update(Ts... vs) { return update(make_numeric_vector({vs...})); }

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
        numeric_matrix<N_, N_> &F_;        // state-transition model
        const numeric_matrix<N_, L_> &B_;  // control-input model
        const numeric_matrix<M_, N_> &H_;  // measurement model
        numeric_matrix<N_, N_> &Q_;        // covariance of the process noise
        numeric_matrix<M_, M_> &R_;        // covariance of the measurement noise
        numeric_vector<N_> x_;             // state vector
        numeric_matrix<N_, N_> P_;         // state covariance, self-initialized as Q_
        const real_t alpha_;               // EMA Smoothing factor for R
        const real_t beta_;                // EMA Smoothing factor for Q

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
         * @param alpha EMA Smoothing factor for R
         * @param beta EMA Smoothing factor for Q
         */
        constexpr adaptive_kalman_filter_t(
                numeric_matrix<N_, N_> &F_matrix,
                const numeric_matrix<N_, L_> &B_matrix,
                const numeric_matrix<M_, N_> &H_matrix,
                numeric_matrix<N_, N_> &Q_matrix,
                numeric_matrix<M_, M_> &R_matrix,
                const numeric_vector<N_> &x_0,
                const real_t &alpha = 0.1,
                const real_t &beta  = 0.1)
            : F_{F_matrix}, B_{B_matrix}, H_{H_matrix},
              Q_{Q_matrix}, R_{R_matrix}, x_{x_0}, P_{Q_matrix},
              alpha_{alpha}, beta_{beta} {}

        constexpr adaptive_kalman_filter_t(const adaptive_kalman_filter_t &) = default;

        constexpr adaptive_kalman_filter_t(adaptive_kalman_filter_t &&) noexcept = default;

        adaptive_kalman_filter_t &operator=(const adaptive_kalman_filter_t &) = default;

        adaptive_kalman_filter_t &operator=(adaptive_kalman_filter_t &&) noexcept = default;

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

        adaptive_kalman_filter_t &operator<<(const numeric_vector<M_> &z) {
            return predict().update(z);
        }

        template<typename... Ts>
        adaptive_kalman_filter_t &update(Ts... vs) { return update(make_numeric_vector({vs...})); }

        const numeric_vector<N_> &state_vector = x_;
        const real_t &state                    = x_[0];
        const numeric_matrix<M_, M_> &R        = R_;
        const numeric_matrix<N_, N_> &Q        = Q_;

    private:
        void adapt_R(const numeric_matrix<M_, M_> &y_yT, const numeric_matrix<M_, M_> &S) {
            R_ = (1 - alpha_) * R_ + alpha_ * (y_yT + S);
        }

        void adapt_Q(const numeric_matrix<N_, M_> &K, const numeric_matrix<M_, M_> &y_yT) {
            Q_ = (1 - beta_) * Q_ + beta_ * (K * y_yT * K.transpose());
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
        using state_func_t           = numeric_vector<N_> (*)(const numeric_vector<N_> &x, const numeric_vector<L_> &u);
        using state_jacobian_t       = numeric_matrix<N_, N_> (*)(const numeric_vector<N_> &x, const numeric_vector<L_> &u);
        using observation_func_t     = numeric_vector<M_> (*)(const numeric_vector<N_> &x);
        using observation_jacobian_t = numeric_matrix<M_, N_> (*)(const numeric_vector<N_> &x);

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

        constexpr extended_kalman_filter_t(const extended_kalman_filter_t &) = default;

        constexpr extended_kalman_filter_t(extended_kalman_filter_t &&) noexcept = default;

        extended_kalman_filter_t &operator=(const extended_kalman_filter_t &) = default;

        extended_kalman_filter_t &operator=(extended_kalman_filter_t &&) noexcept = default;

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

        extended_kalman_filter_t &operator<<(const numeric_vector<M_> &z) {
            return predict().update(z);
        }

        template<typename... Ts>
        extended_kalman_filter_t &update(Ts... vs) { return update(make_numeric_vector({vs...})); }

        const numeric_vector<N_> &state_vector = x_;

        const real_t &state = x_[0];
    };

    template<size_t StateVectorDimension, size_t MeasurementVectorDimension, size_t ControlVectorDimension>
    class adaptive_extended_kalman_filter_t {
    private:
        // Note that numeric_matrix<N_, M_> maps from R_^M_ to R_^N_
        static constexpr size_t N_ = StateVectorDimension;        // ALias
        static constexpr size_t M_ = MeasurementVectorDimension;  // Alias
        static constexpr size_t L_ = ControlVectorDimension;      // Alias

    public:
        using state_func_t           = numeric_vector<N_> (*)(const numeric_vector<N_> &x, const numeric_vector<L_> &u);
        using state_jacobian_t       = numeric_matrix<N_, N_> (*)(const numeric_vector<N_> &x, const numeric_vector<L_> &u);
        using observation_func_t     = numeric_vector<M_> (*)(const numeric_vector<N_> &x);
        using observation_jacobian_t = numeric_matrix<M_, N_> (*)(const numeric_vector<N_> &x);

    protected:
        const state_func_t f_;             // state-transition model
        const state_jacobian_t Fj_;        // state-transition Jacobian
        const observation_func_t h_;       // measurement model
        const observation_jacobian_t Hj_;  // measurement Jacobian
        numeric_matrix<N_, N_> &Q_;        // covariance of the process noise
        numeric_matrix<M_, M_> &R_;        // covariance of the measurement noise
        numeric_vector<N_> x_;             // state vector
        numeric_matrix<N_, N_> P_;         // state covariance, self-initialized as Q_
        const real_t alpha_;               // EMA Smoothing factor for R
        const real_t beta_;                // EMA Smoothing factor for Q

    public:
        constexpr adaptive_extended_kalman_filter_t(
                const state_func_t f_vec_func,
                const state_jacobian_t Fj_mat_func,
                const observation_func_t h_vec_func,
                const observation_jacobian_t Hj_mat_func,
                numeric_matrix<N_, N_> &Q_matrix,
                numeric_matrix<M_, M_> &R_matrix,
                const numeric_vector<N_> &x_0,
                const real_t &alpha = 0.1,
                const real_t &beta  = 0.1)
            : f_(f_vec_func), Fj_{Fj_mat_func}, h_{h_vec_func}, Hj_{Hj_mat_func},
              Q_{Q_matrix}, R_{R_matrix}, x_{x_0}, P_{Q_matrix},
              alpha_{alpha}, beta_{beta} {}

        constexpr adaptive_extended_kalman_filter_t(const adaptive_extended_kalman_filter_t &) = default;

        constexpr adaptive_extended_kalman_filter_t(adaptive_extended_kalman_filter_t &&) noexcept = default;

        adaptive_extended_kalman_filter_t &operator=(const adaptive_extended_kalman_filter_t &) = default;

        adaptive_extended_kalman_filter_t &operator=(adaptive_extended_kalman_filter_t &&) noexcept = default;

        adaptive_extended_kalman_filter_t &predict(const numeric_vector<L_> &u = {}) {
            x_                        = vt::move(f_(x_, u));
            numeric_matrix<N_, N_> F_ = vt::move(Fj_(x_, u));
            P_                        = vt::move(F_ * P_.matmul_T(F_) + Q_);
            return *this;
        }

        adaptive_extended_kalman_filter_t &update(const numeric_vector<M_> &z) {
            numeric_vector<M_> y_          = vt::move(z - h_(x_));
            numeric_matrix<M_, N_> Hjx_    = vt::move(Hj_(x_));
            numeric_matrix<N_, M_> P_Hjx_t = vt::move(P_.matmul_T(Hjx_));
            numeric_matrix<M_, M_> S_      = vt::move(Hjx_ * P_Hjx_t + R_);
            numeric_matrix<N_, M_> K_      = vt::move(P_Hjx_t * S_.inverse());

            x_ += K_ * y_;
            P_ = vt::move((numeric_matrix<N_, N_>::identity() - K_ * Hj_(x_)) * P_);

            numeric_matrix<M_, 1> y_mat = y_.as_matrix_col();
            numeric_matrix<M_, M_> y_yT = y_mat.matmul_T(y_mat);

            adapt_R(y_yT, S_);
            adapt_Q(K_, y_yT);

            return *this;
        }

        adaptive_extended_kalman_filter_t &operator<<(const numeric_vector<M_> &z) {
            return predict().update(z);
        }

        template<typename... Ts>
        adaptive_extended_kalman_filter_t &update(const Ts &...vs) { return update(make_numeric_vector({vs...})); }

        const numeric_vector<N_> &state_vector = x_;
        const real_t &state                    = x_[0];
        const numeric_matrix<M_, M_> &R        = R_;
        const numeric_matrix<N_, N_> &Q        = Q_;

    private:
        void adapt_R(const numeric_matrix<M_, M_> &y_yT, const numeric_matrix<M_, M_> &S) {
            R_ = (1 - alpha_) * R_ + alpha_ * (y_yT + S);
        }

        void adapt_Q(const numeric_matrix<N_, M_> &K, const numeric_matrix<M_, M_> &y_yT) {
            Q_ = (1 - beta_) * Q_ + beta_ * (K * y_yT * K.transpose());
        }
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
    using KF =
            kalman_filter_t<StateVectorDimension, MeasurementVectorDimension, ControlVectorDimension>;

    template<size_t StateVectorDimension, size_t MeasurementVectorDimension, size_t ControlVectorDimension>
    using EKF =
            extended_kalman_filter_t<StateVectorDimension, MeasurementVectorDimension, ControlVectorDimension>;

    namespace basic {
        class KalmanFilter_1D {
        private:
            real_t m_x;  // Estimated state
            real_t m_P;  // Estimated error covariance
            real_t m_Q;  // Process noise covariance
            real_t m_R;  // Measurement noise covariance
            real_t m_K;  // Kalman gain

        public:
            constexpr KalmanFilter_1D() : KalmanFilter_1D(initial_x, initial_P, initial_noise, initial_noise) {}

            constexpr KalmanFilter_1D(const real_t &initial_x, const real_t &initial_P,
                                      const real_t &Q, const real_t &R)
                : m_x(initial_x), m_P(initial_P), m_Q(Q), m_R(R), m_K(0.0) {
            }

            KalmanFilter_1D &predict(const real_t & = 0.0) {
                m_P = m_P + m_Q;
                return *this;
            }

            KalmanFilter_1D &update(const real_t &z) {
                m_K = m_P / (m_P + m_R);
                m_x = m_x + m_K * (z - m_x);
                m_P = (1 - m_K) * m_P;
                return *this;
            }

            KalmanFilter_1D &operator<<(const real_t &z) {
                return predict().update(z);
            }

            [[nodiscard]] constexpr real_t x() const {
                return m_x;
            }

            [[nodiscard]] constexpr real_t P() const {
                return m_P;
            }

            void operator>>(real_t &targ) const {
                targ = x();
            }

            static constexpr real_t initial_x     = 0.0;
            static constexpr real_t initial_P     = 1.0;
            static constexpr real_t initial_noise = 0.1;
        };
    }  // namespace basic
}  // namespace vt

#endif  //VT_LINALG_KALMAN_H
