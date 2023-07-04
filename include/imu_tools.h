#ifndef VT_LINALG_IMU_TOOLS_H
#define VT_LINALG_IMU_TOOLS_H

#include "vt_linalg"

namespace vt {
    typedef struct {
        real_t roll;
        real_t pitch;
        real_t yaw;
    } euler_t;

    typedef struct {
        real_t r;
        real_t i;
        real_t j;
        real_t k;
    } quaternion_t;

    namespace imu_tools {
        euler_t quaternion_to_euler(const quaternion_t &quaternion, bool degrees = true) {
            real_t sqr = quaternion.r * quaternion.r;
            real_t sqi = quaternion.i * quaternion.i;
            real_t sqj = quaternion.j * quaternion.j;
            real_t sqk = quaternion.k * quaternion.k;

            euler_t euler = {
                    atan2(2.0 * (quaternion.i * quaternion.j + quaternion.k * quaternion.r), (sqi - sqj - sqk + sqr)),
                    asin(-2.0 * (quaternion.i * quaternion.k - quaternion.j * quaternion.r) / (sqi + sqj + sqk + sqr)),
                    atan2(2.0 * (quaternion.j * quaternion.k + quaternion.i * quaternion.r), (-sqi - sqj + sqk + sqr)),
            };

            if (degrees) {
                euler.roll *= RAD_TO_DEG;
                euler.pitch *= RAD_TO_DEG;
                euler.yaw *= RAD_TO_DEG;
            }

            return euler;
        }

        numeric_vector<3> quaternion_to_euler(const numeric_vector<4> &quaternion, bool degrees = true) {
            real_t sqr = quaternion[0] * quaternion[0];
            real_t sqi = quaternion[1] * quaternion[1];
            real_t sqj = quaternion[2] * quaternion[2];
            real_t sqk = quaternion[3] * quaternion[3];

            numeric_vector<3> euler;
            euler = {
                    atan2(2.0 * (quaternion[1] * quaternion[2] + quaternion[3] * quaternion[0]),
                          (sqi - sqj - sqk + sqr)),
                    asin(-2.0 * (quaternion[1] * quaternion[3] - quaternion[2] * quaternion[0]) /
                         (sqi + sqj + sqk + sqr)),
                    atan2(2.0 * (quaternion[2] * quaternion[3] + quaternion[1] * quaternion[0]),
                          (-sqi - sqj + sqk + sqr)),
            };

            if (degrees) euler *= RAD_TO_DEG;

            return euler;
        }

        quaternion_t euler_to_quaternion(const euler_t &euler, bool degrees = true) {
            double yaw, pitch, roll;

            // if the angles are in degrees convert them to radians
            if (degrees) {
                yaw = euler.yaw * DEG_TO_RAD;
                pitch = euler.pitch * DEG_TO_RAD;
                roll = euler.roll * DEG_TO_RAD;
            } else {
                yaw = euler.yaw;
                pitch = euler.pitch;
                roll = euler.roll;
            }

            double cy = cos(yaw * 0.5);
            double sy = sin(yaw * 0.5);
            double cp = cos(pitch * 0.5);
            double sp = sin(pitch * 0.5);
            double cr = cos(roll * 0.5);
            double sr = sin(roll * 0.5);

            quaternion_t q;
            q.r = cy * cp * cr + sy * sp * sr;
            q.i = cy * cp * sr - sy * sp * cr;
            q.j = sy * cp * sr + cy * sp * cr;
            q.k = sy * cp * cr - cy * sp * sr;

            return q;
        }

        numeric_vector<4> euler_to_quaternion(const numeric_vector<3> &euler, bool degrees = true) {
            real_t yaw, pitch, roll;

            if (degrees) {
                roll = euler[0] * DEG_TO_RAD;
                pitch = euler[1] * DEG_TO_RAD;
                yaw = euler[2] * DEG_TO_RAD;
            } else {
                yaw = euler[0];
                pitch = euler[1];
                roll = euler[2];
            }

            real_t cy = cos(yaw * 0.5);
            real_t sy = sin(yaw * 0.5);
            real_t cp = cos(pitch * 0.5);
            real_t sp = sin(pitch * 0.5);
            real_t cr = cos(roll * 0.5);
            real_t sr = sin(roll * 0.5);

            numeric_vector<4> q;
            q[0] = cy * cp * cr + sy * sp * sr;
            q[1] = cy * cp * sr - sy * sp * cr;
            q[2] = sy * cp * sr + cy * sp * cr;
            q[3] = sy * cp * cr - cy * sp * sr;

            return q;
        }
    }
}

#endif //VT_LINALG_IMU_TOOLS_H
