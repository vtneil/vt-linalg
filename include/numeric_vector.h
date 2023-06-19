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

    /**
     * Numeric vector template class where the dimension must be known at compile-time
     * and can't be changed by any ways during runtime to prevent unexpected
     * behavior or errors. This avoid heap usages for safety-critical and
     * performance-critical systems.
     *
     * @tparam T data type
     * @tparam Size vector dimension
     */
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
        /**
         * Default constructor, initializes to zero
         */
        constexpr numeric_vector_static_t() = default;

        /**
         * Fill constructor, initializes to fill value
         *
         * @param fill Fill value
         */
        explicit numeric_vector_static_t(T fill) { allocate_fill(fill); }

        /**
         * Copy constructor
         *
         * @param other Other vector
         */
        numeric_vector_static_t(const numeric_vector_static_t &other) { allocate_from(other); }

        /**
         * Move constructor
         *
         * @param other Other vector
         */
        numeric_vector_static_t(numeric_vector_static_t &&other) noexcept {
            for (size_t i = 0; i < Size; ++i) arr_[i] = vt::move(other.arr_[i]);
        }

        /**
         * Array constructor, construct from array
         *
         * @param array Array of entries
         */
        explicit numeric_vector_static_t(const T (&array)[Size]) { allocate_from(array); }

        __VT_FORCE_INLINE T &operator[](size_t index) { return *(arr_ + index); }

        __VT_FORCE_INLINE constexpr const T &operator[](size_t index) const { return *(arr_ + index); }

        __VT_FORCE_INLINE T &at(size_t index) { return operator[](index); };

        __VT_FORCE_INLINE constexpr const T &at(size_t index) const { return operator[](index); };

        __VT_FORCE_INLINE T &operator()(size_t index) { return at(index); }

        __VT_FORCE_INLINE constexpr const T &operator()(size_t index) const { return at(index); }

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

        /**
         * Calculates inner product of this vector (LHS transposed) and the other vector of the same dimension (RHS).
         *
         * @param other Other vector
         * @return Inner product
         */
        T dot(const numeric_vector_static_t &other) const {
            T acc = 0;
            for (size_t i = 0; i < Size; ++i) acc += arr_[i] * other.arr_[i];
            return acc;
        }

        /**
         * Calculates inner product of this vector (LHS transposed) and the other vector of the same dimension (RHS).
         *
         * @param array Other vector as array
         * @return Inner product
         */
        T dot(const T (&array)[Size]) const {
            T acc = 0;
            for (size_t i = 0; i < Size; ++i) acc += arr_[i] * array[i];
            return acc;
        }

        /**
         * Calculates inner product of this vector (LHS transposed) and the other vector of the same dimension (RHS).
         *
         * @param other Other vector
         * @return Inner product
         */
        constexpr T inner(const numeric_vector_static_t &other) const { return dot(other); }

        /**
         * Calculates inner product of this vector (LHS transposed) and the other vector of the same dimension (RHS).
         *
         * @param array Other vector as array
         * @return Inner product
         */
        constexpr T inner(const T (&array)[Size]) const { return dot(array); }

        /**
         * Calculates outer product of this vector (LHS) and the other vector (RHS transposed).
         *
         * @tparam OSize Other vector's dimension
         * @param other Other vector
         * @return Outer product
         */
        template<size_t OSize>
        numeric_matrix_static_t<T, Size, OSize> outer(const numeric_vector_static_t<T, OSize> &other) const {
            numeric_matrix_static_t<T, Size, OSize> result;
            for (size_t i = 0; i < Size; ++i)
                for (size_t j = 0; j < OSize; ++j)
                    result[i][j] = arr_[i] * other.arr_[j];
            return result;
        }

        /**
         * Calculates outer product of this vector (LHS) and the other vector as array (RHS transposed).
         *
         * @tparam OSize Other vector's dimension
         * @param array Other vector as array
         * @return Outer product
         */
        template<size_t OSize>
        numeric_matrix_static_t<T, Size, OSize> outer(const T (&array)[OSize]) const {
            numeric_matrix_static_t<T, Size, OSize> result;
            for (size_t i = 0; i < Size; ++i)
                for (size_t j = 0; j < OSize; ++j)
                    result[i][j] = arr_[i] * array[j];
            return result;
        }

        /**
         * Finds a sum of all entries.
         *
         * @return A sum of all entries
         */
        T sum() const {
            T acc = 0;
            for (size_t i = 0; i < Size; ++i) acc += arr_[i];
            return acc;
        }

        /**
         * Finds a norm of this vector.
         *
         * @return A norm of this vector
         */
        constexpr T norm() const { return pow(dot(*this), 0.5); }

        /**
         * Returns a normalized vector of this vector.
         *
         * @return A normalized vector of this vector
         */
        constexpr numeric_vector_static_t normalize() const { return numeric_vector_static_t(*this) / norm(); }

        /**
         * Checks equality of this vector and the other vector.
         *
         * @tparam OSize Other vector's dimension
         * @param other Other vector to compare
         * @return
         */
        template<size_t OSize>
        bool operator==(const numeric_vector_static_t<T, OSize> &other) const {
            if (this == &other) return true;
            if (Size != OSize) return false;
            for (size_t i = 0; i < Size; ++i) if (arr_[i] != other.arr_[i]) return false;
            return true;
        }

        /**
         * Checks equality of this vector and the other array.
         *
         * @tparam OSize Other array's dimension
         * @param array Array to compare
         * @return
         */
        template<size_t OSize>
        bool operator==(const T (&array)[OSize]) const {
            if (Size != OSize) return false;
            for (size_t i = 0; i < Size; ++i) if (arr_[i] != array[i]) return false;
            return true;
        }

        /**
         * Checks inequality of this vector and the other vector.
         *
         * @tparam OSize Other vector's dimension
         * @param other Other vector to compare
         * @return
         */
        template<size_t OSize>
        constexpr bool operator!=(const numeric_vector_static_t<T, OSize> &other) const { return !operator==(other); }

        /**
         * Checks inequality of this vector and the other array.
         *
         * @tparam OSize Other array's dimension
         * @param array Array to compare
         * @return
         */
        template<size_t OSize>
        constexpr bool operator!=(const T (&array)[OSize]) const { return !operator==(array); }

        /**
         * Checks equality of this vector and the other vector.
         *
         * @tparam OSize Other vector's dimension
         * @param other Other vector to compare
         * @return
         */
        template<size_t OSize>
        constexpr bool equals(const numeric_vector_static_t<T, OSize> &other) const { return operator==(other); }

        /**
         * Checks equality of this vector and the other vector with float/double threshold using
         * equation abs(x - y) < threshold for equality.
         *
         * @tparam OSize Other vector's dimension
         * @param other Other vector to compare
         * @param threshold Equality threshold
         * @return
         */
        template<size_t OSize>
        bool float_equals(const numeric_vector_static_t<T, OSize> &other, real_t threshold = 0.001) const {
            if (this == &other) return true;
            if (Size != OSize) return false;
            for (size_t i = 0; i < Size; ++i) if (abs(arr_[i] - other.arr_[i]) > threshold) return false;
            return true;
        }

        /**
         * Checks equality of this vector and the other array.
         *
         * @tparam OSize Other array's dimension
         * @param array Array to compare
         * @return
         */
        template<size_t OSize>
        constexpr bool equals(const T (&array)[OSize]) const { return operator==(array); }

        /**
         * Returns an iterator to vector's first dimension entry.
         *
         * @return An iterator to vector's first dimension entry
         */
        iterator<T> begin() { return iterator<T>(arr_); }

        /**
         * Returns an iterator to vector's first dimension entry.
         *
         * @return An iterator to vector's first dimension entry
         */
        constexpr iterator<T> begin() const { return iterator<T>(arr_); }

        /**
         * Returns an iterator to vector's last dimension entry.
         *
         * @return An iterator to vector's last dimension entry
         */
        iterator<T> end() { return iterator<T>(arr_ + Size); }

        /**
         * Returns an iterator to vector's last dimension entry.
         *
         * @return An iterator to vector's last dimension entry
         */
        constexpr iterator<T> end() const { return iterator<T>(arr_ + Size); }

        /**
         * Returns vector's dimension (size).
         *
         * @return Vector's dimension (size)
         */
        constexpr size_t size() const { return Size; }

        /**
         * Returns vector's dimension (size).
         *
         * @return Vector's dimension (size)
         */
        constexpr size_t dim() const { return Size; }

        /**
         * Returns vector's dimension (size).
         *
         * @return Vector's dimension (size)
         */
        constexpr size_t dimension() const { return Size; }

        /**
         * Swaps values of this vector with the other vector.
         *
         * @param other Other vector
         */
        void swap(numeric_vector_static_t &other) {
            for (size_t i = 0; i < Size; ++i) vt::swap(arr_[i], other.arr_[i]);
        }

        /**
         * Creates a copy of this vector.
         *
         * @return A copy of this vector
         */
        constexpr numeric_vector_static_t copy() const { return numeric_vector_static_t(*this); }

    private:
        void allocate_zero() { allocate_fill(T()); }

        void allocate_fill(T fill) { vt::fill(arr_, arr_ + Size, fill); }

        void allocate_from(const numeric_vector_static_t &other) {
            static_cast<void>(vt::copy(other.arr_, other.arr_ + Size, arr_));
        }

        void allocate_from(const T (&array)[Size]) {
            static_cast<void>(vt::copy(array, array + Size, arr_));
        }

    public:
        /**
         * Creates a zero matrix.
         *
         * @return Zero matrix
         */
        static constexpr numeric_vector_static_t zeros() { return numeric_vector_static_t(); }

        /**
         * Creates a 1-filled matrix.
         *
         * @return 1-filled matrix
         */
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

    /**
     * Numeric vector where the dimension must be known at compile-time
     * and can't be changed by any ways during runtime to prevent unexpected
     * behavior or errors. This avoid heap usages for safety-critical and
     * performance-critical systems.
     * \n
     * This numeric vector uses real type (real_t), defined as double).
     *
     * @tparam Size Vector dimension
     */
    template<size_t Size>
    using numeric_vector = numeric_vector_static_t<real_t, Size>;

    /**
     * Creates a numeric vector.
     *
     * @tparam Size Vector dimension
     * @param array Array of data
     * @return
     */
    template<size_t Size>
    constexpr numeric_vector<Size> make_numeric_vector(const real_t (&array)[Size]) {
        return numeric_vector<Size>(array);
    }
}

#endif //VNET_LINALG_NUMERIC_VECTOR_H
