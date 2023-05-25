#ifndef VNET_VECTOR_H
#define VNET_VECTOR_H

#include "utils.h"

template<typename T>
class Iterator;

template<typename T>
class Matrix;

template<typename T>
class Vector {
private:
    T *arr_;
    size_t size_;

public:
    Vector() = delete;

    explicit Vector(size_t size) : size_(size) {
        arr_ = new T[size_]();
    }

    Vector(size_t size, T fill) : size_(size) {
        arr_ = new T[size_];
        for (size_t i = 0; i < size_; ++i) arr_[i] = fill;
    }

    Vector(const Vector &other) : size_(other.size_) {
        arr_ = new T[size_];
        memcpy(arr_, other.arr_, other.size_ * sizeof(T));
    }

    Vector(const Iterator<T> &begin, const Iterator<T> &end) : size_(end - begin) {
        arr_ = new T[size_];
        memcpy(arr_, &*begin, size_ * sizeof(T));
    }

    Vector(T *begin, T *end) : size_(end - begin) {
        arr_ = new T[size_];
        memcpy(arr_, begin, size_ * sizeof(T));
    }

    ~Vector() { delete[] arr_; }

    T &operator[](size_t index) { return *(arr_ + index); }

    T &at(size_t index) { return operator[](index); };

    Vector &operator=(const Vector &other) {
        *this = Vector(other);
        return *this;
    }

    Vector &operator+=(const Vector &other) {
        for (size_t i = 0; i < size_; ++i) arr_[i] += other.arr_[i];
        return *this;
    }

    Vector operator+(const Vector &other) const {
        Vector tmp(*this);
        tmp.operator+=(other);
        return tmp;
    }

    Vector add(const Vector &other) const { return operator+(other); }

    Vector &operator-=(const Vector &other) {
        for (size_t i = 0; i < size_; ++i) arr_[i] -= other.arr_[i];
        return *this;
    }

    Vector operator-(const Vector &other) const {
        Vector tmp(*this);
        tmp.operator-=(other);
        return tmp;
    }

    Vector subtract(const Vector &other) const { return operator-(other); }

    T operator*(const Vector &other) const {
        T result = 0;
        for (size_t i = 0; i < size_; ++i) result += arr_[i] * other.arr_[i];
        return result;
    }

    T dot(const Vector &other) const { return operator*(other); }

    T inner(const Vector &other) const { return operator*(other); }

    Matrix<T> outer(const Vector &other) const {

    }

    bool operator==(const Vector &other) const {
        if (size_ != other.size_) return false;
        for (size_t i = 0; i < size_; ++i) if (arr_[i] != other.arr_[i]) return false;
        return true;
    }

    bool operator!=(const Vector &other) const {
        return !operator==(other);
    }

    Iterator<T> begin() { return Iterator<T>(arr_); }

    Iterator<T> end() { return Iterator<T>(arr_ + size_); }

    T &front() { return operator[](0); }

    T &back() { return operator[](size_ - 1); }

    size_t size() const { return size_; }

private:
    void put_array(size_t index, T value) {
        arr_[index] = value;
    }

    template<typename... Ts>
    void put_array(size_t index, T value, Ts ...values) {
        arr_[index] = value;
        put_array(index + 1, values...);
    }

public:
    template<typename... Ts>
    static Vector from(Ts... values) {
        Vector tmp = Vector(sizeof...(values));
        tmp.put_array(0, values...);
        return tmp;
    }
};

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

    bool operator==(const Iterator &other) { return ptr_ == other.ptr_; }

    bool operator!=(const Iterator &other) { return !operator==(other); }

    bool operator==(T *ptr) { return ptr_ == ptr; }

    bool operator!=(T *ptr) { return !operator==(ptr); }
};

#endif //VNET_VECTOR_H
