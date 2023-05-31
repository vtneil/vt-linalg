#ifndef VNET_VECTOR_H
#define VNET_VECTOR_H

#include "utils.h"
//#include <initializer_list>

template<typename T>
class Iterator;

template<typename T>
class Matrix;

template<typename T>
class Vector {
private:
    friend class Matrix<T>;

private:
    T *arr_;
    size_t size_;

public:
    Vector() : size_(0), arr_(nullptr) {}

    explicit Vector(size_t size) : size_(size) { allocate_zero(); }

    Vector(size_t size, T fill) : size_(size) { allocate_fill(fill); }

    Vector(const Vector &other) : size_(other.size_) { allocate_from(other); }

    Vector(Vector &&other) noexcept: size_(other.size_), arr_(other.arr_) {
        other.size_ = 0;
        other.arr_ = nullptr;
    }

    template<size_t N>
    explicit Vector(const T (&array)[N]) : size_(N) { allocate_from(array); }

    Vector(const Iterator<T> &begin, const Iterator<T> &end) : size_(end - begin) {
        arr_ = new T[size_];
        for (Iterator<T> it = begin; it != end; ++it) arr_[it - begin] = *it;
    }

    Vector(T *begin, T *end) : size_(end - begin) {
        arr_ = new T[size_];
        for (size_t i = 0; i < size_; ++i) arr_[i] = *(begin + i);
    }

    ~Vector() { deallocate(); }

    T &operator[](size_t index) { return *(arr_ + index); }

    const T &operator[](size_t index) const { return *(arr_ + index); }

    T &at(size_t index) { return operator[](index); };

    const T &at(size_t index) const { return operator[](index); };

    T &operator()(size_t index) { return at(index); }

    const T &operator()(size_t index) const { return at(index); }

    Vector &operator=(const Vector &other) {
        if (this != &other) {
            deallocate();
            allocate_from(other);
        }
        return *this;
    }

    Vector &operator=(Vector &&other) noexcept {
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
    Vector &operator=(const T (&array)[N]) {
        deallocate();
        allocate_from(array);
        return *this;
    }

    Vector &operator+=(const Vector &other) {
        for (size_t i = 0; i < size_; ++i) arr_[i] += other.arr_[i];
        return *this;
    }

    template<size_t N>
    Vector &operator+=(const T (&array)[N]) {
        for (size_t i = 0; i < size_; ++i) arr_[i] += array[i];
        return *this;
    }

    Vector operator+(const Vector &other) const {
        Vector tmp(*this);
        tmp.operator+=(other);
        return tmp;
    }

    template<size_t N>
    Vector operator+(const T (&array)[N]) const {
        Vector tmp(*this);
        tmp.operator+=(array);
        return tmp;
    }

    Vector add(const Vector &other) const { return operator+(other); }

    template<size_t N>
    Vector add(const T (&array)[N]) const { return operator+(array); }

    Vector &operator-=(const Vector &other) {
        for (size_t i = 0; i < size_; ++i) arr_[i] -= other.arr_[i];
        return *this;
    }

    template<size_t N>
    Vector &operator-=(const T (&array)[N]) {
        for (size_t i = 0; i < size_; ++i) arr_[i] -= array[i];
        return *this;
    }

    Vector operator-(const Vector &other) const {
        Vector tmp(*this);
        tmp.operator-=(other);
        return tmp;
    }

    template<size_t N>
    Vector operator-(const T (&array)[N]) const {
        Vector tmp(*this);
        tmp.operator-=(array);
        return tmp;
    }

    Vector subtract(const Vector &other) const { return operator-(other); }

    template<size_t N>
    Vector subtract(const T (&array)[N]) const { return operator-(array); }

    Vector &operator*=(T rhs) {
        for (size_t i = 0; i < size_; ++i) arr_[i] *= rhs;
        return *this;
    }

    Vector operator*(T rhs) const {
        Vector tmp(*this);
        tmp.operator*=(rhs);
        return tmp;
    }

    Vector &operator/=(T rhs) {
        for (size_t i = 0; i < size_; ++i) arr_[i] /= rhs;
        return *this;
    }

    Vector operator/(T rhs) const {
        Vector tmp(*this);
        tmp.operator/=(rhs);
        return tmp;
    }

    T dot(const Vector &other) const {
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

    constexpr T inner(const Vector &other) const { return dot(other); }

    template<size_t N>
    constexpr T inner(const T (&array)[N]) const { return dot(array); }

    Matrix<T> outer(const Vector &other) const {
        Matrix<T> result(size_, other.size_);
        for (size_t i = 0; i < size_; ++i)
            for (size_t j = 0; j < other.size_; ++j)
                result[i][j] = arr_[i] * other.arr_[j];
        return result;
    }

    template<size_t N>
    Matrix<T> outer(const T (&array)[N]) const {
        Matrix<T> result(size_, N);
        for (size_t i = 0; i < size_; ++i)
            for (size_t j = 0; j < N; ++j)
                result[i][j] = arr_[i] * array[j];
        return result;
    }

    T norm() const { return pow(dot(*this), 0.5); }

    Vector normalize() const { return Vector(*this) / norm(); }

    bool operator==(const Vector &other) const {
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

    bool operator!=(const Vector &other) const { return !operator==(other); }

    template<size_t N>
    bool operator!=(const T (&array)[N]) const { return !operator==(array); }

    bool equals(const Vector &other) const { return operator==(other); }

    bool float_equals(const Vector &other) const {
        if (this == &other) return true;
        if (size_ != other.size_) return false;
        for (size_t i = 0; i < size_; ++i) if (abs(arr_[i] - other.arr_[i]) > 0.001) return false;
        return true;
    }

    template<size_t N>
    bool equals(const T (&array)[N]) const { return operator==(array); }

    void swap(Vector &other) {
        swap_val(arr_, other.arr_);
        swap_val(size_, other.size_);
    }

    Iterator<T> begin() { return Iterator<T>(arr_); }

    Iterator<T> begin() const { return Iterator<T>(arr_); }

    Iterator<T> end() { return Iterator<T>(arr_ + size_); }

    Iterator<T> end() const { return Iterator<T>(arr_ + size_); }

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
        for (size_t i = 0; i < size_; ++i) arr_[i] = fill;
    }

    void allocate_from(const Vector &other) {
        size_ = other.size_;
        arr_ = new T[size_];
        for (size_t i = 0; i < size_; ++i) arr_[i] = other.arr_[i];
    }

    template<size_t N>
    void allocate_from(const T (&array)[N]) {
        size_ = N;
        arr_ = new T[size_];
        for (size_t i = 0; i < size_; ++i) arr_[i] = array[i];
    }

    void deallocate() { delete[] arr_; }

public:
    static Vector zero(size_t n) { return Vector(n); }

    template<typename... Ts>
    static Vector from(Ts... values) {
        Vector tmp = Vector(sizeof...(values));
        tmp.put_array(0, values...);
        return tmp;
    }
};

template<typename T>
Vector<T> operator*(T lhs, const Vector<T> &rhs) {
    Vector<T> tmp(rhs);
    tmp.operator*=(lhs);
    return tmp;
}

template<typename T, size_t N>
Vector<T> operator+(const T (&lhs)[N], const Vector<T> &rhs) {
    Vector<T> tmp(rhs);
    tmp.operator+=(lhs);
    return tmp;
}

template<typename T, size_t N>
Vector<T> operator-(const T (&lhs)[N], const Vector<T> &rhs) {
    Vector<T> tmp(rhs);
    tmp.operator-=(lhs);
    return tmp;
}

template<typename T>
class Iterator {
private:
    T *ptr_;

public:
    explicit Iterator(T *ptr) : ptr_(ptr) {}

    Iterator(const Iterator &other) : ptr_(other.ptr_) {}

    T &operator*() { return *ptr_; }

    T *operator->() { return ptr_; }

    Iterator operator+(size_t n) { return Iterator(ptr_ + n); }

    Iterator &operator++() {
        ++ptr_;
        return *this;
    }

    Iterator operator++(int) {
        Iterator tmp(*this);
        operator++();
        return tmp;
    }

    Iterator operator-(size_t n) { return Iterator(ptr_ - n); }

    size_t operator-(Iterator other) { return ptr_ - other.ptr_; }

    Iterator &operator--() {
        --ptr_;
        return *this;
    }

    Iterator operator--(int) {
        Iterator tmp(*this);
        operator--();
        return tmp;
    }

    bool operator==(const Iterator &other) const { return ptr_ == other.ptr_; }

    bool operator!=(const Iterator &other) const { return !operator==(other); }

    bool operator==(T *ptr) const { return ptr_ == ptr; }

    bool operator!=(T *ptr) const { return !operator==(ptr); }
};

template<typename T>
bool operator==(T *ptr, const Iterator<T> &iterator) { return iterator.operator==(ptr); }

template<typename T>
bool operator!=(T *ptr, const Iterator<T> &iterator) { return !operator==(ptr, iterator); }

using numeric_vector = Vector<double>;

#endif //VNET_VECTOR_H
