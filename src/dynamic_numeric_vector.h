/**
 * @file dynamic_numeric_vector.h
 * @author Vivatsathorn Thitasirivit
 * @date 31 May 2023
 * @brief Numeric vector library
 */

#ifndef VNET_LINALG_DYNAMIC_NUMERIC_VECTOR_H
#define VNET_LINALG_DYNAMIC_NUMERIC_VECTOR_H

#include "standard_utility.h"
#include "iterator.h"

namespace vt {
    template<typename T>
    class numeric_matrix_dynamic_t;

    template<typename T>
    class numeric_vector_dynamic_t {
    private:
        friend class numeric_matrix_dynamic_t<T>;

    private:
        size_t size_;
        T *arr_;

    public:
        numeric_vector_dynamic_t() : size_(0), arr_(nullptr) {}

        explicit numeric_vector_dynamic_t(size_t size) : size_(size) { allocate_zero(); }

        numeric_vector_dynamic_t(size_t size, T fill) : size_(size) { allocate_fill(fill); }

        numeric_vector_dynamic_t(const numeric_vector_dynamic_t &other) : size_(other.size_) { allocate_from(other); }

        numeric_vector_dynamic_t(numeric_vector_dynamic_t &&other) noexcept: size_(other.size_), arr_(other.arr_) {
            other.size_ = 0;
            other.arr_ = nullptr;
        }

        template<size_t N>
        explicit numeric_vector_dynamic_t(const T (&array)[N]) : size_(N) { allocate_from(array); }

        numeric_vector_dynamic_t(const iterator<T> &begin, const iterator<T> &end) : size_(end - begin) {
            arr_ = new T[size_];
            static_cast<void>(vt::copy(begin, end, arr_));
        }

        numeric_vector_dynamic_t(T *begin, T *end) : size_(end - begin) {
            arr_ = new T[size_];
            static_cast<void>(vt::copy(begin, end, arr_));
        }

        ~numeric_vector_dynamic_t() { deallocate(); }

        T &operator[](size_t index) { return *(arr_ + index); }

        const T &operator[](size_t index) const { return *(arr_ + index); }

        T &at(size_t index) { return operator[](index); };

        const T &at(size_t index) const { return operator[](index); };

        T &operator()(size_t index) { return at(index); }

        const T &operator()(size_t index) const { return at(index); }

        numeric_vector_dynamic_t &operator=(const numeric_vector_dynamic_t &other) {
            if (this != &other) {
                deallocate();
                allocate_from(other);
            }
            return *this;
        }

        numeric_vector_dynamic_t &operator=(numeric_vector_dynamic_t &&other) noexcept {
            if (this != &other) {
                deallocate();
                size_ = other.size_;
                arr_ = other.arr_;
                other.size_ = 0;
                other.arr_ = nullptr;
            }
            return *this;
        }

        template<size_t N>
        numeric_vector_dynamic_t &operator=(const T (&array)[N]) {
            deallocate();
            allocate_from(array);
            return *this;
        }

        numeric_vector_dynamic_t &operator+=(const numeric_vector_dynamic_t &other) {
            for (size_t i = 0; i < size_; ++i) arr_[i] += other.arr_[i];
            return *this;
        }

        template<size_t N>
        numeric_vector_dynamic_t &operator+=(const T (&array)[N]) {
            for (size_t i = 0; i < size_; ++i) arr_[i] += array[i];
            return *this;
        }

        numeric_vector_dynamic_t operator+(const numeric_vector_dynamic_t &other) const {
            numeric_vector_dynamic_t tmp(*this);
            tmp.operator+=(other);
            return tmp;
        }

        template<size_t N>
        numeric_vector_dynamic_t operator+(const T (&array)[N]) const {
            numeric_vector_dynamic_t tmp(*this);
            tmp.operator+=(array);
            return tmp;
        }

        numeric_vector_dynamic_t add(const numeric_vector_dynamic_t &other) const { return operator+(other); }

        template<size_t N>
        numeric_vector_dynamic_t add(const T (&array)[N]) const { return operator+(array); }

        numeric_vector_dynamic_t &operator-=(const numeric_vector_dynamic_t &other) {
            for (size_t i = 0; i < size_; ++i) arr_[i] -= other.arr_[i];
            return *this;
        }

        template<size_t N>
        numeric_vector_dynamic_t &operator-=(const T (&array)[N]) {
            for (size_t i = 0; i < size_; ++i) arr_[i] -= array[i];
            return *this;
        }

        numeric_vector_dynamic_t operator-(const numeric_vector_dynamic_t &other) const {
            numeric_vector_dynamic_t tmp(*this);
            tmp.operator-=(other);
            return tmp;
        }

        template<size_t N>
        numeric_vector_dynamic_t operator-(const T (&array)[N]) const {
            numeric_vector_dynamic_t tmp(*this);
            tmp.operator-=(array);
            return tmp;
        }

        numeric_vector_dynamic_t subtract(const numeric_vector_dynamic_t &other) const { return operator-(other); }

        template<size_t N>
        numeric_vector_dynamic_t subtract(const T (&array)[N]) const { return operator-(array); }

        numeric_vector_dynamic_t &operator*=(T rhs) {
            for (size_t i = 0; i < size_; ++i) arr_[i] *= rhs;
            return *this;
        }

        numeric_vector_dynamic_t operator*(T rhs) const {
            numeric_vector_dynamic_t tmp(*this);
            tmp.operator*=(rhs);
            return tmp;
        }

        numeric_vector_dynamic_t &operator/=(T rhs) {
            for (size_t i = 0; i < size_; ++i) arr_[i] /= rhs;
            return *this;
        }

        numeric_vector_dynamic_t operator/(T rhs) const {
            numeric_vector_dynamic_t tmp(*this);
            tmp.operator/=(rhs);
            return tmp;
        }

        T dot(const numeric_vector_dynamic_t &other) const {
            T acc = 0;
            for (size_t i = 0; i < size_; ++i) acc += arr_[i] * other.arr_[i];
            return acc;
        }

        template<size_t N>
        T dot(const T (&array)[N]) const {
            T acc = 0;
            for (size_t i = 0; i < size_; ++i) acc += arr_[i] * array[i];
            return acc;
        }

        constexpr T inner(const numeric_vector_dynamic_t &other) const { return dot(other); }

        template<size_t N>
        constexpr T inner(const T (&array)[N]) const { return dot(array); }

        numeric_matrix_dynamic_t<T> outer(const numeric_vector_dynamic_t &other) const {
            numeric_matrix_dynamic_t<T> result(size_, other.size_);
            for (size_t i = 0; i < size_; ++i)
                for (size_t j = 0; j < other.size_; ++j)
                    result[i][j] = arr_[i] * other.arr_[j];
            return result;
        }

        template<size_t N>
        numeric_matrix_dynamic_t<T> outer(const T (&array)[N]) const {
            numeric_matrix_dynamic_t<T> result(size_, N);
            for (size_t i = 0; i < size_; ++i)
                for (size_t j = 0; j < N; ++j)
                    result[i][j] = arr_[i] * array[j];
            return result;
        }

        T norm() const { return pow(dot(*this), 0.5); }

        numeric_vector_dynamic_t normalize() const { return numeric_vector_dynamic_t(*this) / norm(); }

        bool operator==(const numeric_vector_dynamic_t &other) const {
            if (this == &other) return true;
            if (size_ != other.size_) return false;
            for (size_t i = 0; i < size_; ++i) if (arr_[i] != other.arr_[i]) return false;
            return true;
        }

        template<size_t N>
        bool operator==(const T (&array)[N]) const {
            if (size_ != N) return false;
            for (size_t i = 0; i < size_; ++i) if (arr_[i] != array[i]) return false;
            return true;
        }

        bool operator!=(const numeric_vector_dynamic_t &other) const { return !operator==(other); }

        template<size_t N>
        bool operator!=(const T (&array)[N]) const { return !operator==(array); }

        bool equals(const numeric_vector_dynamic_t &other) const { return operator==(other); }

        bool float_equals(const numeric_vector_dynamic_t &other) const {
            if (this == &other) return true;
            if (size_ != other.size_) return false;
            for (size_t i = 0; i < size_; ++i) if (abs(arr_[i] - other.arr_[i]) > 0.001) return false;
            return true;
        }

        template<size_t N>
        bool equals(const T (&array)[N]) const { return operator==(array); }

        void swap(numeric_vector_dynamic_t &other) {
            vt::swap(arr_, other.arr_);
            vt::swap(size_, other.size_);
        }

        iterator<T> begin() { return iterator<T>(arr_); }

        iterator<T> begin() const { return iterator<T>(arr_); }

        iterator<T> end() { return iterator<T>(arr_ + size_); }

        iterator<T> end() const { return iterator<T>(arr_ + size_); }

        T &front() { return operator[](0); }

        T &back() { return operator[](size_ - 1); }

        size_t size() const { return size_; }

    private:
        void put_array(size_t index, T value) { arr_[index] = value; }

        template<typename... Ts>
        void put_array(size_t index, T value, Ts ...values) {
            arr_[index] = value;
            put_array(index + 1, values...);
        }

        void allocate_zero() { arr_ = new T[size_](); }

        void allocate_fill(T fill) {
            arr_ = new T[size_];
            vt::fill(arr_, arr_ + size_, fill);
        }

        void allocate_from(const numeric_vector_dynamic_t &other) {
            size_ = other.size_;
            arr_ = new T[size_];
            static_cast<void>(vt::copy(other.arr_, other.arr_ + other.size_, arr_));
        }

        template<size_t N>
        void allocate_from(const T (&array)[N]) {
            size_ = N;
            arr_ = new T[size_];
            static_cast<void>(vt::copy(array, array + N, arr_));
        }

        void deallocate() { delete[] arr_; }

    public:
        static numeric_vector_dynamic_t zero(size_t n) { return numeric_vector_dynamic_t(n); }

        template<typename... Ts>
        static numeric_vector_dynamic_t from(Ts... values) {
            numeric_vector_dynamic_t tmp = numeric_vector_dynamic_t(sizeof...(values));
            tmp.put_array(0, values...);
            return tmp;
        }
    };

    template<typename T>
    numeric_vector_dynamic_t<T> operator*(T lhs, const numeric_vector_dynamic_t<T> &rhs) {
        numeric_vector_dynamic_t<T> tmp(rhs);
        tmp.operator*=(lhs);
        return tmp;
    }

    template<typename T, size_t N>
    numeric_vector_dynamic_t<T> operator+(const T (&lhs)[N], const numeric_vector_dynamic_t<T> &rhs) {
        numeric_vector_dynamic_t<T> tmp(rhs);
        tmp.operator+=(lhs);
        return tmp;
    }

    template<typename T, size_t N>
    numeric_vector_dynamic_t<T> operator-(const T (&lhs)[N], const numeric_vector_dynamic_t<T> &rhs) {
        numeric_vector_dynamic_t<T> tmp(rhs);
        tmp.operator-=(lhs);
        return tmp;
    }

    using numeric_vector = numeric_vector_dynamic_t<double>;
}

#endif //VNET_LINALG_DYNAMIC_NUMERIC_VECTOR_H
