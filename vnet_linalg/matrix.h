#ifndef VNET_MATRIX_H
#define VNET_MATRIX_H

#include "utils.h"
#include "vector.h"

template<typename T>
class Matrix {
private:
    T **arr_;
    size_t r_;
    size_t c_;

public:
    Matrix() = delete;

    Matrix(size_t rows, size_t cols) : r_(rows), c_(cols) {
        arr_ = new T *[r_]();
        for (size_t i = 0; i < r_; ++i) arr_[i] = new T[c_]();
    }

    Matrix(size_t rows, size_t cols, T fill) : r_(rows), c_(cols) {
        arr_ = new T *[r_]();
        for (size_t i = 0; i < r_; ++i) {
            arr_[i] = new T[c_]();
            for (size_t j = 0; j < c_; ++j) arr_[i][j] = fill;
        }
    }

    Matrix(const Matrix &other) : r_(other.r_), c_(other.c_) {
        arr_ = new T *[r_]();
        for (size_t i = 0; i < r_; ++i) {
            arr_[i] = new T[c_]();
            for (size_t j = 0; j < c_; ++j) arr_[i][j] = other.arr_[i][j];
        }
    }

    ~Matrix() {
        for (size_t i = 0; i < r_; ++i) delete[] arr_[i];
        delete[] arr_;
    }

    T *&operator[](size_t index) { return *(arr_ + index); }

    T &at(size_t r_index, size_t c_index) { return operator[](r_index)[c_index]; };

    Matrix &operator=(const Matrix &other) {
        *this = Matrix(other);
        return *this;
    }

    Matrix &operator+=(const Matrix &other) {
        for (size_t i = 0; i < r_; ++i) {
            for (size_t j = 0; j < c_; ++j) {
                arr_[i][j] += other.arr_[i][j];
            }
        }
        return *this;
    }

    Matrix operator+(const Matrix &other) const {
        Matrix tmp(*this);
        tmp.operator+=(other);
        return tmp;
    }

    Matrix add(const Matrix &other) const { return operator+(other); }

    Matrix &operator-=(const Matrix &other) {
        for (size_t i = 0; i < r_; ++i) {
            for (size_t j = 0; j < c_; ++j) {
                arr_[i][j] -= other.arr_[i][j];
            }
        }
        return *this;
    }

    Matrix operator-(const Matrix &other) const {
        Matrix tmp(*this);
        tmp.operator-=(other);
        return tmp;
    }

    Matrix subtract(const Matrix &other) const { return operator-(other); }

    Matrix &operator*=(const Matrix &other) {

    }

    Matrix operator*(const Matrix &other) const {

    }

    Matrix matmul(const Matrix &other) const { return operator*(other); }

    bool operator==(const Matrix &other) const {

    }

    bool operator!=(const Matrix &other) const {

    }

    bool can_multiply_with(const Matrix &other) const { return c_ == other.r_; }

    bool can_multiply_with(const Vector<T> &other) const { return c_ == other.size(); }

    Vector<T> row(size_t r_index) const {
        Vector<T> tmp(c_);
        for (size_t i = 0; i < c_; ++i) tmp[i] = arr_[r_index][i];
        return tmp;
    }

    Vector<T> col(size_t c_index) const {
        Vector<T> tmp(r_);
        for (size_t i = 0; i < r_; ++i) tmp[i] = arr_[i][c_index];
        return tmp;
    }

    Vector<T> diag() {
        size_t min_dim = min_val(r_, c_);
        Vector<T> tmp(min_dim);
        for (size_t i = 0; i < min_dim; ++i) tmp[i] = arr_[i][i];
        return tmp;
    }

    Matrix transpose() const {
        Matrix tmp(c_, r_);
    }

    size_t r() const { return r_; }

    size_t c() const { return c_; }

    size_t n() const { return r_ * c_; }
};

#endif //VNET_MATRIX_H
