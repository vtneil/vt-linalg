#ifndef VT_LINALG_COMPLEX_NUMBER_H
#define VT_LINALG_COMPLEX_NUMBER_H

#include "standard_utility.h"

namespace vt {
    class complex_t {
    private:
        real_t real_ = 0.;
        real_t imag_ = 0.;

    public:
        constexpr complex_t() = default;

        constexpr complex_t(real_t real, real_t imag) : real_(real), imag_(imag) {}

        constexpr complex_t(real_t real) : real_(real) {}

        constexpr complex_t(const complex_t &other) = default;

        constexpr complex_t(complex_t &&other) noexcept = default;

        complex_t &operator=(complex_t other) {
            vt::swap(real_, other.real_);
            vt::swap(imag_, other.imag_);
            return *this;
        }

        complex_t &operator=(real_t real) {
            real_ = real;
            imag_ = 0.;
            return *this;
        }

        complex_t &operator+=(const complex_t &rhs) {
            real_ += rhs.real_;
            imag_ += rhs.imag_;
            return *this;
        }

        complex_t &operator+=(real_t rhs) {
            real_ += rhs;
            return *this;
        }

        complex_t operator+(const complex_t &rhs) const {
            complex_t tmp(*this);
            tmp.operator+=(rhs);
            return tmp;
        }

        complex_t operator+(real_t rhs) const {
            complex_t tmp(*this);
            tmp.operator+=(rhs);
            return tmp;
        }

        complex_t &operator++() {
            real_ += 1.;
            return *this;
        }

        const complex_t operator++(int) {
            complex_t tmp(*this);
            operator++();
            return tmp;
        }

        complex_t &operator-=(const complex_t &rhs) {
            real_ -= rhs.real_;
            imag_ -= rhs.imag_;
            return *this;
        }

        complex_t &operator-=(real_t rhs) {
            real_ -= rhs;
            return *this;
        }

        complex_t operator-(const complex_t &rhs) const {
            complex_t tmp(*this);
            tmp.operator-=(rhs);
            return tmp;
        }

        complex_t operator-(real_t rhs) const {
            complex_t tmp(*this);
            tmp.operator-=(rhs);
            return tmp;
        }

        complex_t &operator--() {
            real_ -= 1.;
            return *this;
        }

        const complex_t operator--(int) {
            complex_t tmp(*this);
            tmp.operator--();
            return tmp;
        }

        constexpr complex_t operator-() const { return {-real_, -imag_}; }

        complex_t &operator*=(const complex_t &rhs) {
            real_t temp_real = real_;
            real_            = temp_real * rhs.real_ - imag_ * rhs.imag_;
            imag_            = temp_real * rhs.imag_ + imag_ * rhs.real_;
            return *this;
        }

        complex_t &operator*=(real_t rhs) {
            real_ *= rhs;
            imag_ *= rhs;
            return *this;
        }

        complex_t operator*(const complex_t &rhs) const {
            complex_t tmp(*this);
            tmp.operator*=(rhs);
            return tmp;
        }

        complex_t operator*(real_t rhs) const {
            complex_t tmp(*this);
            tmp.operator*=(rhs);
            return tmp;
        }

        complex_t &operator/=(const complex_t &rhs) {
            real_t denominator = rhs.modulus_pow_2();
            real_t temp_real   = real_;
            real_              = (temp_real * rhs.real_ + imag_ * rhs.imag_) / denominator;
            imag_              = (imag_ * rhs.real_ - temp_real * rhs.imag_) / denominator;
            return *this;
        }

        complex_t &operator/=(real_t rhs) {
            real_ /= rhs;
            imag_ /= rhs;
            return *this;
        }

        complex_t operator/(const complex_t &rhs) const {
            complex_t tmp(*this);
            tmp.operator/=(rhs);
            return tmp;
        }

        complex_t operator/(real_t rhs) const {
            complex_t tmp(*this);
            tmp.operator/=(rhs);
            return tmp;
        }

    public:
        real_t &real() { return real_; }

        real_t &imag() { return imag_; }

        constexpr complex_t conjugate() const { return {real_, -imag_}; }

        constexpr real_t modulus() const { return sqrt(modulus_pow_2()); }

        constexpr real_t modulus_pow_2() const { return real_ * real_ + imag_ * imag_; }
    };

    complex_t operator+(real_t lhs, const complex_t &rhs) { return rhs.operator+(lhs); }

    complex_t operator-(real_t lhs, const complex_t &rhs) { return complex_t(lhs).operator-(rhs); }

    complex_t operator*(real_t lhs, const complex_t &rhs) { return rhs.operator*(lhs); }

    complex_t operator/(real_t lhs, const complex_t &rhs) { return complex_t(lhs).operator/(rhs); }
}// namespace vt

#endif//VT_LINALG_COMPLEX_NUMBER_H
