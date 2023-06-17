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
    class MatrixStatic;

    template<typename T, size_t Size>
    class VectorStatic {
    public:
        static_assert(Size > 0, "Size must be greater than 0.");

    private:
        template<typename U, size_t V, size_t W>
        friend
        class MatrixStatic;

    private:
        T arr_[Size] = {};

    public:
        VectorStatic() = default;

        explicit VectorStatic(T fill) { allocate_fill(fill); }

        VectorStatic(const VectorStatic &other) { allocate_from(other); }

        VectorStatic(VectorStatic &&other) noexcept { for (size_t i = 0; i < Size; ++i) arr_[i] = move(other.arr_[i]); }

        explicit VectorStatic(const T (&array)[Size]) { allocate_from(array); }

        T &operator[](size_t index) { return *(arr_ + index); }

        const T &operator[](size_t index) const { return *(arr_ + index); }

        T &at(size_t index) { return operator[](index); };

        const T &at(size_t index) const { return operator[](index); };

        T &operator()(size_t index) { return at(index); }

        const T &operator()(size_t index) const { return at(index); }

        VectorStatic &operator=(const VectorStatic &other) {
            if (this != &other) allocate_from(other);
            return *this;
        }

        VectorStatic &operator=(VectorStatic &&other) noexcept {
            if (this != &other) for (size_t i = 0; i < Size; ++i) arr_[i] = move(other.arr_[i]);
            return *this;
        }

        VectorStatic &operator=(const T (&array)[Size]) {
            allocate_from(array);
            return *this;
        }

        VectorStatic &operator+=(const VectorStatic &other) {
            for (size_t i = 0; i < Size; ++i) arr_[i] += other.arr_[i];
            return *this;
        }

        VectorStatic &operator+=(const T (&array)[Size]) {
            for (size_t i = 0; i < Size; ++i) arr_[i] += array[i];
            return *this;
        }

        VectorStatic operator+(const VectorStatic &other) const {
            VectorStatic tmp(*this);
            tmp.operator+=(other);
            return tmp;
        }

        VectorStatic operator+(const T (&array)[Size]) const {
            VectorStatic tmp(*this);
            tmp.operator+=(array);
            return tmp;
        }

        VectorStatic add(const VectorStatic &other) const { return operator+(other); }

        VectorStatic add(const T (&array)[Size]) const { return operator+(array); }

        VectorStatic &operator-=(const VectorStatic &other) {
            for (size_t i = 0; i < Size; ++i) arr_[i] -= other.arr_[i];
            return *this;
        }

        VectorStatic &operator-=(const T (&array)[Size]) {
            for (size_t i = 0; i < Size; ++i) arr_[i] -= array[i];
            return *this;
        }

        VectorStatic operator-(const VectorStatic &other) const {
            VectorStatic tmp(*this);
            tmp.operator-=(other);
            return tmp;
        }

        VectorStatic operator-(const T (&array)[Size]) const {
            VectorStatic tmp(*this);
            tmp.operator-=(array);
            return tmp;
        }

        VectorStatic subtract(const VectorStatic &other) const { return operator-(other); }

        VectorStatic subtract(const T (&array)[Size]) const { return operator-(array); }

        VectorStatic &operator*=(T rhs) {
            for (size_t i = 0; i < Size; ++i) arr_[i] *= rhs;
            return *this;
        }

        VectorStatic operator*(T rhs) const {
            VectorStatic tmp(*this);
            tmp.operator*=(rhs);
            return tmp;
        }

        VectorStatic &operator/=(T rhs) {
            for (size_t i = 0; i < Size; ++i) arr_[i] /= rhs;
            return *this;
        }

        VectorStatic operator/(T rhs) const {
            VectorStatic tmp(*this);
            tmp.operator/=(rhs);
            return tmp;
        }

        T dot(const VectorStatic &other) const {
            T acc = 0;
            for (size_t i = 0; i < Size; ++i) acc += arr_[i] * other.arr_[i];
            return acc;
        }

        T dot(const T (&array)[Size]) const {
            T acc = 0;
            for (size_t i = 0; i < Size; ++i) acc += arr_[i] * array[i];
            return acc;
        }

        T inner(const VectorStatic &other) const { return dot(other); }

        T inner(const T (&array)[Size]) const { return dot(array); }

        template<size_t OSize>
        MatrixStatic<T, Size, OSize> outer(const VectorStatic<T, OSize> &other) const {
            MatrixStatic<T, Size, OSize> result;
            for (size_t i = 0; i < Size; ++i)
                for (size_t j = 0; j < OSize; ++j)
                    result[i][j] = arr_[i] * other.arr_[j];
            return result;
        }

        template<size_t OSize>
        MatrixStatic<T, Size, OSize> outer(const T (&array)[OSize]) const {
            MatrixStatic<T, Size, OSize> result;
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

        T norm() const { return pow(dot(*this), 0.5); }

        VectorStatic normalize() const { return VectorStatic(*this) / norm(); }

        template<size_t OSize>
        bool operator==(const VectorStatic<T, OSize> &other) const {
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
        bool operator!=(const VectorStatic<T, OSize> &other) const { return !operator==(other); }

        template<size_t OSize>
        bool operator!=(const T (&array)[OSize]) const { return !operator==(array); }

        template<size_t OSize>
        bool equals(const VectorStatic<T, OSize> &other) const { return operator==(other); }

        template<size_t OSize>
        bool float_equals(const VectorStatic<T, OSize> &other) const {
            if (this == &other) return true;
            if (Size != OSize) return false;
            for (size_t i = 0; i < Size; ++i) if (abs(arr_[i] - other.arr_[i]) > 0.001) return false;
            return true;
        }

        template<size_t OSize>
        bool equals(const T (&array)[OSize]) const { return operator==(array); }

        void swap(VectorStatic &other) { for (size_t i = 0; i < Size; ++i) vt::swap(arr_[i], other.arr_[i]); }

        Iterator<T> begin() { return Iterator<T>(arr_); }

        Iterator<T> begin() const { return Iterator<T>(arr_); }

        Iterator<T> end() { return Iterator<T>(arr_ + Size); }

        Iterator<T> end() const { return Iterator<T>(arr_ + Size); }

        T &front() { return operator[](0); }

        T &back() { return operator[](Size - 1); }

        constexpr size_t size() const { return Size; }

    private:
        void put_array(size_t index, T value) { arr_[index] = value; }

        template<typename... Ts>
        void put_array(size_t index, T value, Ts ...values) {
            arr_[index] = value;
            if (index < Size - 1) put_array(index + 1, values...);
        }

        void allocate_zero() { allocate_fill(T()); }

        void allocate_fill(T fill) { vt::fill(arr_, arr_ + Size, fill); }

        void allocate_from(const VectorStatic &other) {
            static_cast<void>(vt::copy(other.arr_, other.arr_ + Size, arr_));
        }

        void allocate_from(const T (&array)[Size]) {
            static_cast<void>(vt::copy(array, array + Size, arr_));
        }

    public:
        static VectorStatic zero() { return VectorStatic(); }

        template<typename... Ts>
        static VectorStatic from(Ts... values) {
            VectorStatic tmp = VectorStatic(sizeof...(values));
            tmp.put_array(0, values...);
            return tmp;
        }
    };

    template<typename T, size_t Size>
    VectorStatic<T, Size> operator*(T lhs, const VectorStatic<T, Size> &rhs) {
        VectorStatic<T, Size> tmp(rhs);
        tmp.operator*=(lhs);
        return tmp;
    }

    template<typename T, size_t Size>
    VectorStatic<T, Size> operator+(const T (&lhs)[Size], const VectorStatic<T, Size> &rhs) {
        VectorStatic<T, Size> tmp(rhs);
        tmp.operator+=(lhs);
        return tmp;
    }

    template<typename T, size_t Size>
    VectorStatic<T, Size> operator-(const T (&lhs)[Size], const VectorStatic<T, Size> &rhs) {
        VectorStatic<T, Size> tmp(rhs);
        tmp.operator-=(lhs);
        return tmp;
    }

    template<size_t Size>
    using numeric_vector = vt::VectorStatic<double, Size>;
}

#endif //VNET_LINALG_NUMERIC_VECTOR_H
