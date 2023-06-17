/**
 * @file numeric_vector.h
 * @author Vivatsathorn Thitasirivit
 * @date 15 June 2023
 * @brief Numeric static vector library
 */

#ifndef VNET_LINALG_NUMERIC_VECTOR_H
#define VNET_LINALG_NUMERIC_VECTOR_H

#include "standard_utility.h"
#include "iterator.h"

namespace vt {
    template<typename T, size_t Row, size_t Col>
    class numeric_matrix_static_t;

    template<typename T, size_t Size>
    class numeric_vector_static_t {
    public:
        static_assert(Size > 0, "Capacity must be greater than 0.");

    private:
        template<typename U, size_t V, size_t W>
        friend
        class numeric_matrix_static_t;

    private:
        T arr_[Size] = {};

    public:
        numeric_vector_static_t() = default;

        explicit numeric_vector_static_t(T fill) { allocate_fill(fill); }

        numeric_vector_static_t(const numeric_vector_static_t &other) { allocate_from(other); }

        numeric_vector_static_t(numeric_vector_static_t &&other) noexcept {
            for (size_t i = 0; i < Size; ++i) arr_[i] = vt::move(other.arr_[i]);
        }

        explicit numeric_vector_static_t(const T (&array)[Size]) { allocate_from(array); }

        T &operator[](size_t index) { return *(arr_ + index); }

        constexpr const T &operator[](size_t index) const { return *(arr_ + index); }

        T &at(size_t index) { return operator[](index); };

        constexpr const T &at(size_t index) const { return operator[](index); };

        T &operator()(size_t index) { return at(index); }

        constexpr const T &operator()(size_t index) const { return at(index); }

        numeric_vector_static_t &operator=(const numeric_vector_static_t &other) {
            if (this != &other) allocate_from(other);
            return *this;
        }

        numeric_vector_static_t &operator=(numeric_vector_static_t &&other) noexcept {
            if (this != &other) for (size_t i = 0; i < Size; ++i) arr_[i] = vt::move(other.arr_[i]);
            return *this;
        }

        numeric_vector_static_t &operator=(const T (&array)[Size]) {
            allocate_from(array);
            return *this;
        }

        numeric_vector_static_t &operator+=(const numeric_vector_static_t &other) {
            for (size_t i = 0; i < Size; ++i) arr_[i] += other.arr_[i];
            return *this;
        }

        numeric_vector_static_t &operator+=(const T (&array)[Size]) {
            for (size_t i = 0; i < Size; ++i) arr_[i] += array[i];
            return *this;
        }

        numeric_vector_static_t operator+(const numeric_vector_static_t &other) const {
            numeric_vector_static_t tmp(*this);
            tmp.operator+=(other);
            return tmp;
        }

        numeric_vector_static_t operator+(const T (&array)[Size]) const {
            numeric_vector_static_t tmp(*this);
            tmp.operator+=(array);
            return tmp;
        }

        constexpr numeric_vector_static_t add(const numeric_vector_static_t &other) const { return operator+(other); }

        constexpr numeric_vector_static_t add(const T (&array)[Size]) const { return operator+(array); }

        numeric_vector_static_t &operator-=(const numeric_vector_static_t &other) {
            for (size_t i = 0; i < Size; ++i) arr_[i] -= other.arr_[i];
            return *this;
        }

        numeric_vector_static_t &operator-=(const T (&array)[Size]) {
            for (size_t i = 0; i < Size; ++i) arr_[i] -= array[i];
            return *this;
        }

        numeric_vector_static_t operator-(const numeric_vector_static_t &other) const {
            numeric_vector_static_t tmp(*this);
            tmp.operator-=(other);
            return tmp;
        }

        numeric_vector_static_t operator-(const T (&array)[Size]) const {
            numeric_vector_static_t tmp(*this);
            tmp.operator-=(array);
            return tmp;
        }

        constexpr numeric_vector_static_t subtract(const numeric_vector_static_t &other) const {
            return operator-(other);
        }

        constexpr numeric_vector_static_t subtract(const T (&array)[Size]) const { return operator-(array); }

        numeric_vector_static_t &operator*=(T rhs) {
            for (size_t i = 0; i < Size; ++i) arr_[i] *= rhs;
            return *this;
        }

        numeric_vector_static_t operator*(T rhs) const {
            numeric_vector_static_t tmp(*this);
            tmp.operator*=(rhs);
            return tmp;
        }

        numeric_vector_static_t &operator/=(T rhs) {
            for (size_t i = 0; i < Size; ++i) arr_[i] /= rhs;
            return *this;
        }

        numeric_vector_static_t operator/(T rhs) const {
            numeric_vector_static_t tmp(*this);
            tmp.operator/=(rhs);
            return tmp;
        }

        T dot(const numeric_vector_static_t &other) const {
            T acc = 0;
            for (size_t i = 0; i < Size; ++i) acc += arr_[i] * other.arr_[i];
            return acc;
        }

        T dot(const T (&array)[Size]) const {
            T acc = 0;
            for (size_t i = 0; i < Size; ++i) acc += arr_[i] * array[i];
            return acc;
        }

        constexpr T inner(const numeric_vector_static_t &other) const { return dot(other); }

        constexpr T inner(const T (&array)[Size]) const { return dot(array); }

        template<size_t OSize>
        numeric_matrix_static_t<T, Size, OSize> outer(const numeric_vector_static_t<T, OSize> &other) const {
            numeric_matrix_static_t<T, Size, OSize> result;
            for (size_t i = 0; i < Size; ++i)
                for (size_t j = 0; j < OSize; ++j)
                    result[i][j] = arr_[i] * other.arr_[j];
            return result;
        }

        template<size_t OSize>
        numeric_matrix_static_t<T, Size, OSize> outer(const T (&array)[OSize]) const {
            numeric_matrix_static_t<T, Size, OSize> result;
            for (size_t i = 0; i < Size; ++i)
                for (size_t j = 0; j < OSize; ++j)
                    result[i][j] = arr_[i] * array[j];
            return result;
        }

        T sum() const {
            T acc = 0;
            for (size_t i = 0; i < Size; ++i) acc += arr_[i];
            return acc;
        }

        constexpr T norm() const { return pow(dot(*this), 0.5); }

        constexpr numeric_vector_static_t normalize() const { return numeric_vector_static_t(*this) / norm(); }

        template<size_t OSize>
        bool operator==(const numeric_vector_static_t<T, OSize> &other) const {
            if (this == &other) return true;
            if (Size != OSize) return false;
            for (size_t i = 0; i < Size; ++i) if (arr_[i] != other.arr_[i]) return false;
            return true;
        }

        template<size_t OSize>
        bool operator==(const T (&array)[OSize]) const {
            if (Size != OSize) return false;
            for (size_t i = 0; i < Size; ++i) if (arr_[i] != array[i]) return false;
            return true;
        }

        template<size_t OSize>
        constexpr bool operator!=(const numeric_vector_static_t<T, OSize> &other) const { return !operator==(other); }

        template<size_t OSize>
        constexpr bool operator!=(const T (&array)[OSize]) const { return !operator==(array); }

        template<size_t OSize>
        constexpr bool equals(const numeric_vector_static_t<T, OSize> &other) const { return operator==(other); }

        template<size_t OSize>
        bool float_equals(const numeric_vector_static_t<T, OSize> &other) const {
            if (this == &other) return true;
            if (Size != OSize) return false;
            for (size_t i = 0; i < Size; ++i) if (abs(arr_[i] - other.arr_[i]) > 0.001) return false;
            return true;
        }

        template<size_t OSize>
        constexpr bool equals(const T (&array)[OSize]) const { return operator==(array); }

        iterator<T> begin() { return iterator<T>(arr_); }

        constexpr iterator<T> begin() const { return iterator<T>(arr_); }

        iterator<T> end() { return iterator<T>(arr_ + Size); }

        constexpr iterator<T> end() const { return iterator<T>(arr_ + Size); }

        T &front() { return operator[](0); }

        T &back() { return operator[](Size - 1); }

        constexpr size_t size() const { return Size; }

        void swap(numeric_vector_static_t &other) {
            for (size_t i = 0; i < Size; ++i) vt::swap(arr_[i], other.arr_[i]);
        }

        constexpr numeric_vector_static_t copy() const { return numeric_vector_static_t(*this); }

    private:
        void put_array(size_t index, T value) { arr_[index] = value; }

        template<typename... Ts>
        void put_array(size_t index, T value, Ts ...values) {
            arr_[index] = value;
            if (index < Size - 1) put_array(index + 1, values...);
        }

        void allocate_zero() { allocate_fill(T()); }

        void allocate_fill(T fill) { vt::fill(arr_, arr_ + Size, fill); }

        void allocate_from(const numeric_vector_static_t &other) {
            static_cast<void>(vt::copy(other.arr_, other.arr_ + Size, arr_));
        }

        void allocate_from(const T (&array)[Size]) {
            static_cast<void>(vt::copy(array, array + Size, arr_));
        }

    public:
        static constexpr numeric_vector_static_t zeros() { return numeric_vector_static_t(); }

        static constexpr numeric_vector_static_t ones() { return numeric_vector_static_t(1); }
    };

    template<typename T, size_t Size>
    numeric_vector_static_t<T, Size> operator*(T lhs, const numeric_vector_static_t<T, Size> &rhs) {
        numeric_vector_static_t<T, Size> tmp(rhs);
        tmp.operator*=(lhs);
        return tmp;
    }

    template<typename T, size_t Size>
    numeric_vector_static_t<T, Size> operator+(const T (&lhs)[Size], const numeric_vector_static_t<T, Size> &rhs) {
        numeric_vector_static_t<T, Size> tmp(rhs);
        tmp.operator+=(lhs);
        return tmp;
    }

    template<typename T, size_t Size>
    numeric_vector_static_t<T, Size> operator-(const T (&lhs)[Size], const numeric_vector_static_t<T, Size> &rhs) {
        numeric_vector_static_t<T, Size> tmp(rhs);
        tmp.operator-=(lhs);
        return tmp;
    }

    template<size_t Size>
    using numeric_vector = numeric_vector_static_t<real_t, Size>;

    template<size_t Size>
    constexpr numeric_vector<Size> make_numeric_vector(const real_t (&array)[Size]) {
        return numeric_vector<Size>(array);
    }
}

#endif //VNET_LINALG_NUMERIC_VECTOR_H
