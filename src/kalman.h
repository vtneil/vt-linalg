/**
 * @file kalman.h
 * @author Vivatsathorn Thitasirivit
 * @date 31 May 2023
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

#include "utils.h"
#include "dynamic_numeric_vector.h"
#include "dynamic_numeric_matrix.h"

class SimpleKalmanFilter {
protected:
    double dt_;
    numeric_matrix A_;
    numeric_matrix B_;
    numeric_matrix C_;
    numeric_matrix D_;
    numeric_matrix G_;
    numeric_matrix Q_;
    numeric_matrix R_;

protected:
    numeric_vector x_;
    numeric_vector x_p;
    numeric_vector y_;
    numeric_vector y_p;
    numeric_matrix P_;
    numeric_matrix P_p;
    numeric_matrix S_;
    numeric_matrix K_;

public:
    SimpleKalmanFilter(double dt, numeric_matrix A, numeric_matrix B, numeric_matrix C,
                       numeric_matrix D, numeric_matrix G, numeric_matrix Q, numeric_matrix R) :
            dt_(dt), A_(move(A)), B_(move(B)), C_(move(C)),
            D_(move(D)), G_(move(G)), Q_(move(Q)), R_(move(R)) {

    }

    SimpleKalmanFilter(double dt, numeric_matrix &&A, numeric_matrix &&B, numeric_matrix &&C,
                       numeric_matrix &&D, numeric_matrix &&G, numeric_matrix &&Q, numeric_matrix &&R) :
            dt_(dt), A_(move(A)), B_(move(B)), C_(move(C)),
            D_(move(D)), G_(move(G)), Q_(move(Q)), R_(move(R)) {

    }

    void predict(const numeric_vector &u) {
        x_p = move(A_ * x_ + B_ * u);
        P_p = move(A_ * P_ * A_.transpose() + G_ * Q_ * G_.transpose());
    }

    void update(const numeric_vector &y, const numeric_vector &u) {
        y_p = move(C_ * x_p + D_ * u);
        y_ = move(y - y_p);
        S_ = move(C_ * P_p * C_.transpose() + R_);
        K_ = move(P_p * C_.transpose() * S_.inv());
        x_ = move(x_p + K_ * y_);
        numeric_matrix KC = move(K_ * C_);
        P_ = move((numeric_matrix::id(min_val(KC.r(), KC.c())) - KC) * P_p);
    }
};

/**
 * Alias for Simple Kalman Filter
 */
using KF = SimpleKalmanFilter;
using KalmanFilter = SimpleKalmanFilter;

#endif //VNET_LINALG_KALMAN_H
