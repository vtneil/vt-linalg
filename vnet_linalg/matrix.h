#ifndef VNET_MATRIX_H
#define VNET_MATRIX_H

#include "utils.h"
#include "vector.h"

template<typename T>
class Matrix {
private:
    friend class Vector<T>;

private:
    Vector<Vector<T> *> *vector_;
    size_t r_;
    size_t c_;

public:
    Matrix() = delete;

    explicit Matrix(size_t n) : r_(n), c_(n) {
        allocate_zero();
    }

    Matrix(size_t rows, size_t cols) : r_(rows), c_(cols) {
        allocate_zero();
    }

    Matrix(size_t rows, size_t cols, T fill) : r_(rows), c_(cols) {
        allocate_fill(fill);
    }

    Matrix(const Matrix &other) : r_(other.r_), c_(other.c_) {
        allocate_from(other);
    }

    ~Matrix() {
        deallocate();
    }

    Vector<T> &operator[](size_t index) { return *(vector_->operator[](index)); }

    const Vector<T> &operator[](size_t index) const { return *(vector_->operator[](index)); }

    T &at(size_t r_index, size_t c_index) { return operator[](r_index).operator[](c_index); };

    const T &at(size_t r_index, size_t c_index) const { return operator[](r_index).operator[](c_index); };

    T &operator()(size_t r_index, size_t c_index) { return at(r_index, c_index); }

    const T &operator()(size_t r_index, size_t c_index) const { return at(r_index, c_index); }

    Matrix &operator=(const Matrix &other) {
        if (this != &other) {
            deallocate();
            allocate_from(other);
        }
        return *this;
    }

    Matrix &operator+=(const Matrix &other) {
        for (size_t i = 0; i < r_; ++i) operator[](i) += other[i];
        return *this;
    }

    Matrix operator+(const Matrix &other) const {
        Matrix tmp(*this);
        tmp.operator+=(other);
        return tmp;
    }

    Matrix add(const Matrix &other) const { return operator+(other); }

    Matrix &operator-=(const Matrix &other) {
        for (size_t i = 0; i < r_; ++i) operator[](i) -= other[i];
        return *this;
    }

    Matrix operator-(const Matrix &other) const {
        Matrix tmp(*this);
        tmp.operator-=(other);
        return tmp;
    }

    Matrix subtract(const Matrix &other) const { return operator-(other); }

    Matrix &operator*=(const Matrix &other) {
        Matrix tmp = operator*(other);
        deallocate();
        allocate_from(tmp);
        return *this;
    }

    Matrix operator*(const Matrix &other) const {
        if (max_val(max_val(r_, other.r_), max_val(c_, other.c_)) > 128)
            return mm_strassen(*this, other);
        else
            return mm_naive(*this, other);
    }

    Matrix &operator*=(T rhs) {
        for (auto &x: *vector_) *x *= rhs;
        return *this;
    }

    Matrix operator*(T rhs) const {
        Matrix tmp(*this);
        tmp.operator*=(rhs);
        return tmp;
    }

    /**
     *
     * @param other Vector to be left-multiplied with this matrix
     * @return Multiplied vector b = Ax
     */
    Vector<T> operator*(const Vector<T> &other) const {
        Vector<T> tmp(r_);
        for (size_t i = 0; i < r_; ++i) tmp[i] = *(vector_->operator[](i)) * other;
        return tmp;
    }

    Matrix matmul(const Matrix &other) const { return operator*(other); }

    Matrix &operator^=(size_t n) {
        Matrix tmp = operator^(n);
        deallocate();
        allocate_from(tmp);
        return *this;
    }

    Matrix operator^(size_t n) {
        Matrix product(*this);

        if (n == 1) return product;

        product = product.operator^(n / 2);

        if (n % 2 == 0) return product.operator*(product);
        else return product.operator*(product).operator*(*this);
    }

    Matrix matpow(size_t n) {
        return operator^(n);
    }

    bool operator==(const Matrix &other) const {
        if (this == &other) return true;
        if (r_ != other.r_ || c_ != other.c_) return false;
        for (size_t i = 0; i < r_; ++i)
            for (size_t j = 0; j < c_; ++j)
                if (vector_->operator[](i)->operator[](j) != other.vector_->operator[](i)->operator[](j))
                    return false;
        return true;
    }

    bool operator!=(const Matrix &other) const { return !operator==(other); }

    bool is_square() const { return r_ == c_; }

    bool can_multiply_with(const Matrix &other) const { return c_ == other.r_; }

    bool can_multiply_with(const Vector<T> &other) const { return c_ == other.size_; }

    Vector<T> row(size_t r_index) const {
        return Vector<T>(operator[](r_index));
    }

    Vector<T> col(size_t c_index) const {
        Vector<T> tmp(r_);
        for (size_t i = 0; i < r_; ++i) tmp[i] = vector_[i][c_index];
        return tmp;
    }

    Vector<T> diag() {
        size_t min_dim = min_val(r_, c_);
        Vector<T> tmp(min_dim);
        for (size_t i = 0; i < min_dim; ++i) tmp[i] = vector_[i][i];
        return tmp;
    }

    Matrix transpose() const {
        Matrix tmp(c_, r_);
        for (size_t i = 0; i < r_; ++i) {
            for (size_t j = 0; j < c_; ++j) {
                tmp.vector_->operator[](j)->operator[](i) = vector_->operator[](i)->operator[](j);
            }
        }
        return tmp;
    }

    T det() const {
        Pair<Matrix<T>, Matrix<T>> lu = LU();

        T det = 1;
        size_t r_swaps = 0;
        size_t n = min_val(r_, c_);

        for (size_t i = 0; i < n; ++i) {
            det *= lu.upper()[i][i];
            if (lu.upper()[i][i] == 0) return 0;
            if (lu.lower()[i][i] != 1) ++r_swaps;
        }

        return r_swaps % 2 == 0 ? det : -det;
    }

    T tr() const {
        T acc = 0;
        size_t n = min_val(r_, c_);
        for (size_t i = 0; i < n; ++i) acc += vector_->operator[](i)->operator[](i);
        return acc;
    }

    Matrix<T> inv() const {
        size_t n = min_val(r_, c_);
        Pair<Matrix, Matrix> lu = LU();
        Matrix inverse(n);
        for (size_t i = 0; i < n; ++i) {
            Vector<T> b(n);
            b[i] = 1;
            Vector<T> y = lu_forward_sub(lu.lower(), b);
            Vector<T> x = lu_backward_sub(lu.upper(), y);
            for (size_t j = 0; j < n; ++j) inverse[j][i] = x[j];
        }
        return inverse;
    }

    Matrix<T> RRE() const {

    }

    Pair<Matrix, Matrix> LU() const {
        size_t n = min_val(r_, c_);
        Matrix lower(n, n);
        Matrix upper(n, n);

        for (size_t i = 0; i < n; ++i) {
            for (size_t k = 0; k < n; ++k) {
                T sum_ = 0;
                for (size_t j = 0; j < i; ++j) sum_ += lower[i][j] * upper[j][k];
                upper[i][k] = vector_->operator[](i)->operator[](k) - sum_;
            }
            for (size_t k = i; k < n; ++k) {
                if (i == k) lower[i][i] = 1;
                else {
                    T sum_ = 0;
                    for (size_t j = 0; j < i; ++j) sum_ += lower[k][j] * upper[j][i];
                    lower[k][i] = (vector_->operator[](k)->operator[](i) - sum_) / upper[i][i];
                }
            }
        }
        return {lower, upper};
    }

    void append_row() {

    }

    void append_col() {

    }

    Iterator<Vector<T> *> begin() { return vector_->begin(); }

    Iterator<Vector<T> *> begin() const { return vector_->begin(); }

    Iterator<Vector<T> *> end() { return vector_->end(); }

    Iterator<Vector<T> *> end() const { return vector_->end(); }

    size_t r() const { return r_; }

    size_t c() const { return c_; }

    size_t n() const { return r_ * c_; }

private:
    static Vector<T> lu_forward_sub(const Matrix<double> &L, const Vector<double> &b) {
        size_t n = min_val(L.r_, L.c_);
        Vector<T> y(n);
        for (size_t i = 0; i < n; ++i) {
            y[i] = b[i];
            for (size_t j = 0; j < i; ++j) y[i] -= L[i][j] * y[j];
            y[i] /= L[i][i];
        }
        return y;
    }

    static Vector<T> lu_backward_sub(const Matrix<double> &U, const Vector<double> &y) {
        size_t n = min_val(U.r_, U.c_);
        Vector<T> x(n);
        for (size_t i = n - 1;; --i) {
            x[i] = y[i];
            for (size_t j = i + 1; j < i; ++j) x[i] -= U[i][j] * x[j];
            x[i] /= U[i][i];
            if (i == 0) break;
        }
        return x;
    }

    static Matrix mm_naive(const Matrix &A, const Matrix &B) {
        Matrix C(A.r_, B.c_);
        for (int i = 0; i < A.r_; i++)
            for (int j = 0; j < B.c_; j++)
                for (int k = 0; k < A.c_; k++)
                    C[i][j] += A[i][k] * B[k][j];
        return C;
    }

    static Matrix mm_strassen(const Matrix &A, const Matrix &B) {
        return mm_naive(A, B);
    }

    void put_array(size_t index, Vector<T> &vector) {
        delete vector_->operator[](index);
        vector_->operator[](index) = new Vector<T>(vector);
        c_ = vector.size();
    }

    template<typename... Vectors>
    void put_array(size_t index, Vector<T> &vector, Vectors &...vectors) {
        delete vector_->operator[](index);
        vector_->operator[](index) = new Vector<T>(vector);
        put_array(index + 1, vectors...);
    }

    void allocate_zero() {
        vector_ = new Vector<Vector<T> *>(r_);
        for (size_t i = 0; i < r_; ++i) vector_->operator[](i) = new Vector<T>(c_);
    }

    void allocate_fill(T fill) {
        vector_ = new Vector<Vector<T> *>(r_);
        for (size_t i = 0; i < r_; ++i) {
            vector_->operator[](i) = new Vector<T>(c_);
            for (size_t j = 0; j < c_; ++j) vector_->operator[](i)->operator[](j) = fill;
        }
    }

    void allocate_from(const Matrix &other) {
        r_ = other.r_;
        c_ = other.c_;
        vector_ = new Vector<Vector<T> *>(r_);
        for (size_t i = 0; i < r_; ++i) {
            vector_->operator[](i) = new Vector<T>(c_);
            for (size_t j = 0; j < c_; ++j)
                vector_->operator[](i)->operator[](j) = other.vector_->operator[](i)->operator[](j);
        }
    }

    void deallocate() {
        for (size_t i = 0; i < r_; ++i) delete vector_->operator[](i);
        delete vector_;
    }

public:
    static Matrix zero(size_t r, size_t c) { return Matrix(r, c); }

    static Matrix zero(size_t n) { return Matrix(n); }

    static Matrix identity(size_t n) { return Matrix::diagonal(n, 1); }

    static Matrix I1() { return identity(1); }

    static Matrix I2() { return identity(2); }

    static Matrix I3() { return identity(3); }

    static Matrix I4() { return identity(4); }

    static Matrix I5() { return identity(5); }

    static Matrix diagonal(size_t n, T value) {
        Matrix tmp(n);
        for (size_t i = 0; i < n; ++i) tmp[i][i] = value;
        return tmp;
    }

    template<typename... Vectors>
    static Matrix from(Vectors... vectors) {
        Matrix tmp(sizeof...(vectors), 0);
        tmp.put_array(0, vectors...);
        return tmp;
    }

    static Matrix from_row(const Vector<T> &vector) { return Matrix::from(vector); }

    static Matrix from_col(const Vector<T> &vector) {
        Matrix result(vector.size_, 1);
        for (size_t i = 0; i < vector.size_; ++i) result[i][0] = vector.arr_[i];
        return result;
    }
};

template<typename T>
Matrix<T> operator*(T lhs, const Matrix<T> &rhs) {
    Matrix<T> tmp(rhs);
    tmp.operator*=(lhs);
    return tmp;
}

template<typename T>
T det(const Matrix<T> &A) {
    return A.det();
}

#endif //VNET_MATRIX_H
