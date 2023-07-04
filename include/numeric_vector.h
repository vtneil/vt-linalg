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
    namespace detail {
        template<typename T, size_t Size, size_t IMax, size_t I>
        void helper_assign_to(const T (&array)[Size], T &var) {
            var = array[I];
        }

        template<typename T, size_t Size, size_t IMax, size_t I, typename... Ts>
        void helper_assign_to(const T (&array)[Size], T &var, Ts &...vars) {
            helper_assign_to<T, Size, IMax, I>(array, var);
            if (I < IMax) helper_assign_to<T, Size, IMax, I + 1>(array, vars...);
        }
    }

    namespace impl {
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

        private:
            template<size_t... I>
            constexpr numeric_vector_static_t(const T &fill, vt::index_sequence<I...>)
                    : arr_{(static_cast<void>(I), fill)...} {}

        public:
            /**
             * Fill constructor, initializes to fill value
             *
             * @param fill Fill value
             */
            constexpr explicit numeric_vector_static_t(const T &fill)
                    : numeric_vector_static_t(fill, vt::make_index_sequence<Size>()) {}

            /**
             * Copy constructor
             *
             * @param other Other vector
             */
            constexpr numeric_vector_static_t(const numeric_vector_static_t &other) = default;

            /**
             * Move constructor
             *
             * @param other Other vector
             */
            constexpr numeric_vector_static_t(numeric_vector_static_t &&other) noexcept = default;

        private:
            template<size_t... I>
            constexpr numeric_vector_static_t(const T (&array)[Size], vt::index_sequence<I...>)
                    : arr_{array[I]...} {}

        public:
            /**
             * Array constructor, construct from array
             *
             * @param array Array of entries
             */
            constexpr explicit numeric_vector_static_t(const T (&array)[Size])
                    : numeric_vector_static_t(array, vt::make_index_sequence<Size>()) {}

            /**
             * Extended constructor
             *
             * @tparam S1
             * @tparam S2
             * @param v1
             * @param v2
             */
            template<size_t S1, size_t S2>
            numeric_vector_static_t(const numeric_vector_static_t<T, S1> &v1,
                                    const numeric_vector_static_t<T, S2> &v2) {
                insert<0>(v1);
                insert<S1>(v2);
            }

            /**
             * Extended constructor
             *
             * @tparam S1
             * @tparam S2
             * @param a1
             * @param a2
             */
            template<size_t S1, size_t S2>
            numeric_vector_static_t(const T (&a1)[S1],
                                    const T (&a2)[S2]) {
                insert<0>(a1);
                insert<S1>(a2);
            }

            FORCE_INLINE T &operator[](size_t index) { return *(arr_ + index); }

            FORCE_INLINE constexpr const T &operator[](size_t index) const { return *(arr_ + index); }

            FORCE_INLINE T &at(size_t index) { return operator[](index); };

            FORCE_INLINE constexpr const T &at(size_t index) const { return operator[](index); };

            FORCE_INLINE T &operator()(size_t index) { return at(index); }

            FORCE_INLINE constexpr const T &operator()(size_t index) const { return at(index); }

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

            constexpr numeric_vector_static_t add(const numeric_vector_static_t &other) const {
                return operator+(other);
            }

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
            constexpr bool operator!=(const numeric_vector_static_t<T, OSize> &other) const {
                return !operator==(other);
            }

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
             * equation abs(x_ - y) < threshold for equality.
             *
             * @tparam OSize Other vector's dimension
             * @param other Other vector to compare
             * @param threshold Equality threshold
             * @return
             */
            template<size_t OSize>
            bool float_equals(const numeric_vector_static_t<T, OSize> &other, real_t threshold = 1e-10) const {
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
             * Insert the vector into the vector.
             *
             * @tparam pos Position to insert
             * @tparam OSize
             * @param v Vector to insert
             * @return Reference to this vector
             */
            template<size_t pos = 0, size_t OSize>
            numeric_vector_static_t &insert(const numeric_vector_static_t<T, OSize> &v) {
                static_assert(pos < Size, "Insertion failed! Position must be within range.");
                static_assert(pos + OSize <= Size, "Insertion failed! Vector out of range.");
                for (size_t i = 0; i < OSize; ++i) arr_[pos + i] = v[i];
                return *this;
            }

            /**
             * Insert the vector into the vector.
             *
             * @tparam pos Position to insert
             * @tparam OSize
             * @param array Vector to insert
             * @return Reference to this vector
             */
            template<size_t pos = 0, size_t OSize>
            numeric_vector_static_t &insert(const T (&array)[OSize]) {
                static_assert(pos < Size, "Insertion failed! Position must be within range.");
                static_assert(pos + OSize <= Size, "Insertion failed! Vector out of range.");
                for (size_t i = 0; i < OSize; ++i) arr_[pos + i] = array[i];
                return *this;
            }

            /**
             * Gets a slice of this vector.
             *
             * @tparam from From index
             * @tparam to To index
             * @return Sliced vector
             */
            template<size_t from, size_t to>
            numeric_vector_static_t<T, to - from> slice() {
                static_assert(from < to, "from must be less than to.");
                static_assert(to <= Size, "Slice range is out of range.");
                numeric_vector_static_t<T, to - from> result;
                for (size_t i = 0; i < to - from; ++i) result[i] = arr_[from + i];
                return result;
            }

            /**
             * Converts this vector to a matrix representation of column vector
             *
             * @return Column vector as matrix
             */
            numeric_matrix_static_t<T, Size, 1> as_matrix_col() {
                numeric_matrix_static_t<T, Size, 1> result;
                for (size_t i = 0; i < Size; ++i) result[i][0] = arr_[i];
                return result;
            }

            /**
             * Converts this vector to a matrix representation of row vector
             *
             * @return Row vector as matrix
             */
            numeric_matrix_static_t<T, 1, Size> as_matrix_row() {
                numeric_matrix_static_t<T, 1, Size> result;
                for (size_t i = 0; i < Size; ++i) result[0][i] = arr_[i];
                return result;
            }

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
             * Swaps entries with the other vector.
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

            template<typename... Ts>
            void assign_to(Ts &...vars) {
                static_assert(sizeof...(Ts) <= Size,
                              "Number of variables must not be greater than vector\'s dimension.");
                detail::helper_assign_to<T, Size, sizeof...(Ts), 0>(arr_, vars...);
            }

        private:
            void allocate_zero() { allocate_fill(T()); }

            void allocate_fill(const T &fill) { vt::fill(arr_, arr_ + Size, fill); }

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
    using numeric_vector = impl::numeric_vector_static_t<real_t, Size>;

    /**
     * Helper class for make_numeric_vector function. Not intended for user's usages.
     */
    namespace detail {
        template<typename... Args>
        struct size_sum;

        template<template<typename, size_t> class Class, typename T, size_t Size>
        struct size_sum<Class<T, Size>> {
            static constexpr size_t value = Size;
        };

        template<template<typename, size_t> class Class, typename T, size_t Size, size_t... Sizes>
        struct size_sum<Class<T, Size>, Class<T, Sizes>...> {
            static constexpr size_t value = Size + size_sum<impl::numeric_vector_static_t < T, Sizes>...>::value;
        };
    }

    /**
     * Creates a numeric vector.
     *
     * @tparam Size Vector dimension
     * @param array Array of data
     * @return Numeric vector
     */
    template<size_t Size>
    constexpr numeric_vector<Size> make_numeric_vector(const real_t (&array)[Size]) {
        return numeric_vector<Size>(array);
    }

    /**
     * Creates a numeric vector (copy).
     *
     * @tparam Size Vector dimension
     * @param vector Vector to copy from
     * @return Numeric vector
     */
    template<size_t Size>
    constexpr numeric_vector<Size> make_numeric_vector(const numeric_vector<Size> &vector) {
        return numeric_vector<Size>(vector);
    }

    /**
     * Creates a numeric vector by extending multiple numeric vectors.
     *
     * @tparam S1
     * @tparam S2
     * @param v1 Vector 1
     * @param v2 Vector 2
     * @return Numeric vector
     */
    template<size_t S1, size_t S2>
    constexpr numeric_vector<S1 + S2>
    make_numeric_vector(const numeric_vector<S1> &v1,
                        const numeric_vector<S2> &v2) {
        return numeric_vector<S1 + S2>(v1, v2);
    }

    /**
     * Creates a numeric vector by extending multiple numeric vectors.
     *
     * @tparam S1
     * @tparam S2
     * @tparam Ss
     * @param v1 Vector 1
     * @param v2 Vector 2
     * @param vs Vectors
     * @return Numeric vector
     */
    template<size_t S1, size_t S2, size_t... Ss>
    constexpr numeric_vector<vt::detail::size_sum<numeric_vector<S1>, numeric_vector<S2>, numeric_vector<Ss>...>::value>
    make_numeric_vector(const numeric_vector<S1> &v1,
                        const numeric_vector<S2> &v2,
                        const numeric_vector<Ss> &... vs) {
        return make_numeric_vector(make_numeric_vector(v1, v2), vs...);
    }

    namespace detail {
        using namespace impl;

        template<typename T, size_t Size>
        constexpr numeric_vector_static_t <T, Size> make_numeric_vector_static_t() {
            return numeric_vector_static_t<T, Size>();
        }

        template<typename T, size_t Size>
        constexpr numeric_vector_static_t <T, Size> make_numeric_vector_static_t(const T (&array)[Size]) {
            return numeric_vector_static_t<T, Size>(array);
        }

        template<typename T, size_t S1, size_t S2>
        constexpr numeric_vector_static_t<T, S1 + S2>
        make_numeric_vector_static_t(const T (&a1)[S1],
                                     const T (&a2)[S2]) {
            return numeric_vector_static_t<T, S1 + S2>(a1, a2);
        }

        template<typename T, size_t S1, size_t S2, size_t... Ss>
        constexpr numeric_vector_static_t<T, size_sum<
                numeric_vector_static_t <
                T, S1>, numeric_vector_static_t < T, S2>, numeric_vector_static_t <T, Ss>...>::value>

        make_numeric_vector_static_t(const T (&a1)[S1],
                                     const T (&a2)[S1],
                                     const T (&...as)[Ss]) {
            return make_numeric_vector_static_t(make_numeric_vector_static_t(a1, a2), as...);
        }

        template<typename T>
        constexpr numeric_vector_static_t<T, 1> make_numeric_vector_static_t(const T &val) {
            return numeric_vector_static_t<T, 1>(val);
        }

        template<typename T, typename... Ts>
        constexpr numeric_vector_static_t<T, 1 + sizeof...(Ts)>
        make_numeric_vector_static_t(const T &val, const Ts &...vals) {
            return make_numeric_vector_static_t({val, vals...});
        }

        template<typename T, size_t Row, size_t Col, size_t ...I>
        constexpr numeric_vector_static_t <numeric_vector_static_t<T, Col>, Row>
        make_nested(const T (&array)[Row][Col], vt::index_sequence<I...>) {
            return make_numeric_vector_static_t(make_numeric_vector_static_t(array[I])...);
        }

        template<typename T, size_t Row, size_t Col>
        constexpr numeric_vector_static_t <numeric_vector_static_t<T, Col>, Row>
        make_nested(const T (&array)[Row][Col]) {
            return make_nested(array, vt::make_index_sequence<Row>());
        }
    }
}

#endif //VNET_LINALG_NUMERIC_VECTOR_H
