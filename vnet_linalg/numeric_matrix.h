#ifndef VNET_MATRIX_H
#define VNET_MATRIX_H

#include "utils.h"
#include "numeric_vector.h"

template<typename T_>
class MatrixLU;

template<typename T_>
class Matrix {
private:
    friend class Vector<T_>;

private:
    static constexpr size_t STRASSEN_DIMENSION = 1024;
    static constexpr size_t STRASSEN_THRESHOLD = STRASSEN_DIMENSION * STRASSEN_DIMENSION;
    Vector<Vector<T_>> vector_;
    size_t r_;
    size_t c_;

public:
    Matrix() : r_(0), c_(0), vector_(move(Vector<Vector<T_>>())) {}

    explicit Matrix(size_t n) : r_(n), c_(n) { allocate_zero(); }

    Matrix(size_t rows, size_t cols) : r_(rows), c_(cols) { allocate_zero(); }

    Matrix(size_t rows, size_t cols, T_ fill) : r_(rows), c_(cols) { allocate_fill(fill); }

    Matrix(const Matrix &other) : r_(other.r_), c_(other.c_) { allocate_from(other); }

    Matrix(Matrix &&other) noexcept: r_(other.r_), c_(other.c_) {
        steal(move(other));
        other.r_ = 0;
        other.c_ = 0;
    }

    template<size_t R, size_t C>
    explicit Matrix(const T_ (&array)[R][C]) : r_(R), c_(C) { allocate_from(array); }

    template<size_t R>
    explicit Matrix(const Vector<T_> (&array)[R]) : r_(R), c_(array[0].size_) { allocate_from(array); }

    Matrix(const Matrix &M11, const Matrix &M12, const Matrix &M21, const Matrix &M22) {
        if (M11.r_ == M12.r_ && M21.r_ == M22.r_ && M11.c_ == M21.c_ && M12.c_ == M22.c_) {
            r_ = M11.r_ + M21.r_;
            c_ = M11.c_ + M12.c_;
            allocate_zero();
            insert(0, 0, M11);
            insert(0, c_ / 2, M12);
            insert(r_ / 2, 0, M21);
            insert(r_ / 2, c_ / 2, M22);
        } else {
            r_ = max_val(M11.r_ + M21.r_, M12.r_ + M22.r_);
            c_ = max_val(M11.c_ + M12.c_, M21.c_ + M22.c_);
            allocate_zero();
        }
    }

    Vector<T_> &operator[](size_t index) { return vector_[index]; }

    const Vector<T_> &operator[](size_t index) const { return vector_[index]; }

    T_ &at(size_t r_index, size_t c_index) { return vector_[r_index][c_index]; };

    const T_ &at(size_t r_index, size_t c_index) const { return vector_[r_index][c_index]; };

    T_ &operator()(size_t r_index, size_t c_index) { return at(r_index, c_index); }

    const T_ &operator()(size_t r_index, size_t c_index) const { return at(r_index, c_index); }

    Matrix &operator=(const Matrix &other) {
        if (this != &other) {
            if (r_ == other.r_ && c_ == other.c_) {
                insert(other);
            } else {
                allocate_from(other);
            }
        }
        return *this;
    }

    Matrix &operator=(Matrix &&other) noexcept {
        if (this != &other) {
            steal(move(other));
            other.r_ = 0;
            other.c_ = 0;
        }
        return *this;
    }

    template<size_t R, size_t C>
    Matrix &operator=(const T_ (&array)[R][C]) {
        if (r_ == R && c_ == C) {
            insert(array);
        } else {
            allocate_from(array);
        }
        return *this;
    }

    template<size_t R>
    Matrix &operator=(const Vector<T_> (&vectors)[R]) {
        allocate_from(vectors);
        return *this;
    }

    Matrix &operator+=(const Matrix &other) { return iadd(*this, *this, other); }

    Matrix operator+(const Matrix &other) const {
        Matrix tmp(*this);
        tmp.operator+=(other);
        return tmp;
    }

    Matrix add(const Matrix &other) const { return operator+(other); }

    static Matrix &iadd(Matrix &C, const Matrix &A, const Matrix &B) {
        for (int i = 0; i < C.r_; ++i)
            for (int j = 0; j < C.c_; ++j)
                C[i][j] = A[i][j] + B[i][j];
        return C;
    }

    Matrix &operator-=(const Matrix &other) { return isub(*this, *this, other); }

    Matrix operator-(const Matrix &other) const {
        Matrix tmp(*this);
        tmp.operator-=(other);
        return tmp;
    }

    Matrix sub(const Matrix &other) const { return operator-(other); }

    static Matrix &isub(Matrix &C, const Matrix &A, const Matrix &B) {
        for (int i = 0; i < C.r_; ++i)
            for (int j = 0; j < C.c_; ++j)
                C[i][j] = A[i][j] - B[i][j];
        return C;
    }

    Matrix &operator*=(const Matrix &other) {
        steal(move(operator*(other)));
        return *this;
    }

    Matrix operator*(const Matrix &other) const {
        Matrix C(r_, other.c_);
        return imatmul(C, *this, other);
    }

    Matrix &operator*=(T_ rhs) {
        for (auto &x: vector_) x *= rhs;
        return *this;
    }

    Matrix operator*(T_ rhs) const {
        Matrix tmp(*this);
        tmp.operator*=(rhs);
        return tmp;
    }

    Vector<T_> operator*(const Vector<T_> &other) const {
        Vector<T_> tmp(r_);
        for (size_t i = 0; i < r_; ++i) tmp[i] = vector_[i] * other;
        return tmp;
    }

    Matrix matmul(const Matrix &other) const { return operator*(other); }

    static Matrix &imatmul(Matrix &C, const Matrix &A, const Matrix &B) {
        if (A.r() != 1 &&
            B.r() != 1 &&
            A.c() != 1 &&
            B.c() != 1 &&
            max_val(A.r() * A.c(), B.r() * B.c()) >= STRASSEN_THRESHOLD)
            return mm_strassen(C, A, B);
        return mm_naive(C, A, B);
    }

    Matrix matmul_naive(const Matrix &other) const {
        Matrix C(r_, other.c_);
        return mm_naive(C, *this, other);
    }

    Vector<T_> transform(const Vector<T_> &other) const { return operator*(other); }

    Matrix &operator^=(size_t n) {
        steal(move(operator^(n)));
        return *this;
    }

    Matrix operator^(size_t n) {
        Matrix product(*this);
        if (n == 0) return id(min_val(r_, c_));
        if (n == 1) return product;
        product = move(product.operator^(n / 2));
        if (n % 2 == 0) return product.operator*(product);
        else return product.operator*(product).operator*(*this);
    }

    Matrix matpow(size_t n) { return operator^(n); }

    Matrix matpow_naive(size_t n) {
        Matrix product = move(id(r_));
        for (size_t i = 0; i < n; ++i) product *= *this;
        return product;
    }

    bool operator==(const Matrix &other) const {
        if (this == &other) return true;
        if (r_ != other.r_ || c_ != other.c_) return false;
        for (size_t i = 0; i < r_; ++i)
            if (vector_[i] != other.vector_[i])
                return false;
        return true;
    }

    template<size_t R, size_t C>
    bool operator==(const T_ (&array)[R][C]) const {
        if (r_ != R || c_ != C) return false;
        for (size_t i = 0; i < r_; ++i)
            if (vector_[i] != array[i])
                return false;
        return true;
    }

    bool operator!=(const Matrix &other) const { return !operator==(other); }

    template<size_t R, size_t C>
    bool operator!=(const T_ (&array)[R][C]) const { return !operator==(array); }

    bool equals(const Matrix &other) const { return operator==(other); }

    template<size_t R, size_t C>
    bool equals(const T_ (&array)[R][C]) const { return operator==(array); }

    bool is_square() const { return r_ == c_; }

    bool can_multiply_with(const Matrix &other) const { return c_ == other.r_; }

    bool can_multiply_with(const Vector<T_> &other) const { return c_ == other.size_; }

    Vector<T_> row(size_t r_index) const { return Vector<T_>(operator[](r_index)); }

    Vector<T_> col(size_t c_index) const {
        Vector<T_> tmp(r_);
        for (size_t i = 0; i < r_; ++i) tmp[i] = vector_[i][c_index];
        return tmp;
    }

    Vector<T_> diag() {
        size_t min_dim = min_val(r_, c_);
        Vector<T_> tmp(min_dim);
        for (size_t i = 0; i < min_dim; ++i) tmp[i] = vector_[i][i];
        return tmp;
    }

    Matrix T() const { return transpose(); }

    Matrix transpose() const {
        Matrix tmp(c_, r_);
        for (size_t i = 0; i < r_; ++i)
            for (size_t j = 0; j < c_; ++j)
                tmp.vector_[j][i] = vector_[i][j];
        return tmp;
    }

    T_ det() const {
        MatrixLU<T_> lu = move(LU());
        T_ det = 1;
        size_t r_swaps = 0;
        size_t n = min_val(r_, c_);
        for (size_t i = 0; i < n; ++i) {
            det *= lu.u()[i][i];
            if (lu.u()[i][i] == 0) return 0;
            if (lu.l()[i][i] != 1) ++r_swaps;
        }
        return r_swaps % 2 == 0 ? det : -det;
    }

    T_ tr() const {
        T_ acc = 0;
        size_t n = min_val(r_, c_);
        for (size_t i = 0; i < n; ++i) acc += vector_[i][i];
        return acc;
    }

    Matrix inv() const {
        MatrixLU<T_> lu = move(LU());
        return inv_ut(lu.u()) * inv_lt(lu.l());
    }

    Matrix RRE() const {
        Matrix m(*this);
        size_t lead = 0;
        for (size_t r = 0; r < r_; ++r) {
            if (lead >= c_) return m;
            size_t i;
            for (i = r; m[i][lead] == 0;) {
                ++i;
                if (i == r_) {
                    i = r;
                    ++lead;
                    if (lead == c_) return m;
                }
            }
            m[i].swap(m[r]);
            T_ val = m[r][lead];
            for (size_t j = 0; j < c_; ++j) m[r][j] /= val;
            for (i = 0; i < r_; ++i) {
                if (i != r) {
                    val = m[i][lead];
                    for (size_t j = 0; j < c_; ++j)
                        m[i][j] -= val * m[r][j];
                }
            }
            ++lead;
        }
        m.fix_zero();
        return m;
    }

    MatrixLU<T_> LU() const {
        size_t n = min_val(r_, c_);
        Matrix lower(n, n);
        Matrix upper(n, n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t k = 0; k < n; ++k) {
                T_ sum_ = 0;
                for (size_t j = 0; j < i; ++j) sum_ += lower[i][j] * upper[j][k];
                upper[i][k] = vector_[i][k] - sum_;
            }
            for (size_t k = i; k < n; ++k) {
                if (i == k) lower[i][i] = 1;
                else {
                    T_ sum_ = 0;
                    for (size_t j = 0; j < i; ++j) sum_ += lower[k][j] * upper[j][i];
                    lower[k][i] = (vector_[k][i] - sum_) / upper[i][i];
                }
            }
        }
        return {lower, upper};
    }

    void insert(const Matrix &M) { insert(0, 0, M); }

    template<size_t R, size_t C>
    void insert(const T_ (&array)[R][C]) { insert(0, 0, array); }

    void insert(size_t pos_row, size_t pos_col, const Matrix &M) {
        if (pos_row >= r_ || pos_col >= c_ || pos_row + M.r_ > r_ || pos_col + M.c_ > c_) return;
        for (size_t i = 0; i < M.r_; ++i)
            for (size_t j = 0; j < M.c_; ++j)
                vector_[i + pos_row][j + pos_col] = M[i][j];
    }

    template<size_t R, size_t C>
    void insert(size_t pos_row, size_t pos_col, const T_ (&array)[R][C]) {
        if (pos_row >= r_ || pos_col >= c_ || pos_row + R > r_ || pos_col + C > c_) return;
        for (size_t i = 0; i < R; ++i)
            for (size_t j = 0; j < C; ++j)
                vector_[i + pos_row][j + pos_col] = array[i][j];
    }

    Matrix slice(size_t r1, size_t c1, size_t r2, size_t c2) const {
        if (r1 > r2) swap_val(r1, r2);
        if (c1 > c2) swap_val(c1, c2);
        r2 = min_val(r2, r_);
        c2 = min_val(c2, c_);
        size_t r = r2 - r1;
        size_t c = c2 - c1;
        Matrix result(r, c);
        if (r == 0 || c == 0 || r1 >= r_ || c1 >= c_) return result;
        for (size_t i = 0; i < r; ++i)
            for (size_t j = 0; j < c; ++j)
                result[i][j] = vector_[r1 + i][c1 + j];
        return result;
    }

    Iterator<Vector<T_>> begin() { return vector_.begin(); }

    Iterator<Vector<T_>> begin() const { return vector_.begin(); }

    Iterator<Vector<T_>> end() { return vector_.end(); }

    Iterator<Vector<T_>> end() const { return vector_.end(); }

    size_t r() const { return r_; }

    size_t c() const { return c_; }

    size_t n() const { return r_ * c_; }

public:
    constexpr static bool is_multiple_of_2(size_t x) { return x != 0 && !(x & (x - 1)); }

    static size_t closest_2(size_t x) {
        if (is_multiple_of_2(x)) return x;
        size_t result = 1;
        while (x != 0) {
            x >>= 1;
            result <<= 1;
        }
        return result;
    }

private:
    static Matrix &inv_lt(Matrix &L) {
        size_t n = min_val(L.r_, L.c_);
        Matrix X(n);
        for (size_t k = 0; k < n; ++k) {
            X[k][k] = 1 / L[k][k];
            for (size_t i = k + 1; i < n; ++i) {
                T_ acc = 0;
                for (size_t j = k; j < i; ++j)
                    acc -= L[i][j] * X[j][k];
                X[i][k] = acc / L[i][i];
            }
        }
        L = move(X);
        return L;
    }

    static Matrix inv_ut(const Matrix &U) {
        Matrix tmp = move(U.transpose());
        return inv_lt(tmp).transpose();
    }

    static Matrix &mm_naive(Matrix &C, const Matrix &A, const Matrix &B) {
        for (size_t i = 0; i < A.r_; ++i) {
            Vector<T_> &row_C = C[i];
            for (size_t k = 0; k < A.c_; ++k) {
                const T_ &val_A = A[i][k];
                for (size_t j = 0; j < B.c_; ++j) {
                    row_C[j] += val_A * B[k][j];
                }
            }
        }
        return C;
    }

    static Matrix &mm_strassen(Matrix &C, const Matrix &A, const Matrix &B) {
        size_t rA = closest_2(A.r_);
        size_t cA = closest_2(A.c_);
        size_t rB = closest_2(B.r_);
        size_t cB = closest_2(B.c_);

        Matrix A11, A12, A21, A22, B11, B12, B21, B22;

        if (rA == A.r_ && cA == A.c_ && rB == B.r_ && cB == B.c_) {
            A11 = move(A.slice(0, 0, rA / 2, cA / 2));
            A12 = move(A.slice(0, cA / 2, rA / 2, cA));
            A21 = move(A.slice(rA / 2, 0, rA, cA / 2));
            A22 = move(A.slice(rA / 2, cA / 2, rA, cA));
            B11 = move(B.slice(0, 0, rB / 2, cB / 2));
            B12 = move(B.slice(0, cB / 2, rB / 2, cB));
            B21 = move(B.slice(rB / 2, 0, rB, cB / 2));
            B22 = move(B.slice(rB / 2, cB / 2, rB, cB));
        } else {
            Matrix Ap(rA, cA);
            Ap.insert(A);
            Matrix Bp(rB, cB);
            Bp.insert(B);
            A11 = move(Ap.slice(0, 0, rA / 2, cA / 2));
            A12 = move(Ap.slice(0, cA / 2, rA / 2, cA));
            A21 = move(Ap.slice(rA / 2, 0, rA, cA / 2));
            A22 = move(Ap.slice(rA / 2, cA / 2, rA, cA));
            B11 = move(Bp.slice(0, 0, rB / 2, cB / 2));
            B12 = move(Bp.slice(0, cB / 2, rB / 2, cB));
            B21 = move(Bp.slice(rB / 2, 0, rB, cB / 2));
            B22 = move(Bp.slice(rB / 2, cB / 2, rB, cB));
            C = move(Matrix(rA, cB));
        }

        Matrix M1 = move((A11 + A22) * (B11 + B22));
        Matrix M2 = move((A21 + A22) * B11);
        Matrix M3 = move(A11 * (B12 - B22));
        Matrix M4 = move(A22 * (B21 - B11));
        Matrix M5 = move((A11 + A12) * B22);
        Matrix M6 = move((A21 - A11) * (B11 + B12));
        Matrix M7 = move((A12 - A22) * (B21 + B22));
        C.insert(0, 0, M1 + M4 - M5 + M7);
        C.insert(0, C.c_ / 2, M3 + M5);
        C.insert(C.r_ / 2, 0, M2 + M4);
        C.insert(C.r_ / 2, C.c_ / 2, M1 - M2 + M3 + M6);

        if (rA != A.r_ || cB != B.c_)
            C = C.slice(0, 0, A.r_, B.c_);
        return C;
    }

    void fix_zero() {
        for (size_t i = 0; i < r_; ++i)
            for (size_t j = 0; j < c_; ++j)
                if (vector_[i][j] == 0.0)
                    vector_[i][j] = 0.0;
    }

    void put_array(size_t index, const Vector<T_> &vector) {
        vector_[index] = vector;
        c_ = vector.size();
    }

    template<size_t N>
    void put_array(size_t index, const T_ (&array)[N]) {
        vector_[index] = array;
        c_ = N;
    }

    template<typename... Vectors>
    void put_array(size_t index, const Vector<T_> &vector, const Vectors &...vectors) {
        vector_[index] = vector;
        put_array(index + 1, vectors...);
    }

    template<size_t N, typename... Arrays>
    void put_array(size_t index, const T_ (&array)[N], const Arrays (&...arrays)[N]) {
        vector_[index] = array;
        put_array(index + 1, arrays...);
    }

    void allocate_zero() { vector_ = move(Vector<Vector<T_>>(r_, Vector<T_>(c_))); }

    void allocate_fill(T_ fill) { vector_ = move(Vector<Vector<T_>>(r_, Vector<T_>(c_, fill))); }

    void steal(Matrix &&other) {
        r_ = other.r_;
        c_ = other.c_;
        vector_ = move(other.vector_);
    }

    void allocate_from(const Matrix &other) {
        r_ = other.r_;
        c_ = other.c_;
        vector_ = other.vector_;
    }

    template<size_t R, size_t C>
    void allocate_from(const T_ (&array)[R][C]) {
        r_ = R;
        c_ = C;
        allocate_zero();
        for (size_t i = 0; i < r_; ++i)
            for (size_t j = 0; j < c_; ++j)
                vector_[i][j] = array[i][j];
    }

    template<size_t R>
    void allocate_from(const Vector<T_> (&vectors)[R]) {
        r_ = R;
        c_ = vectors[0].size_;
        allocate_zero();
        for (size_t i = 0; i < r_; ++i)
            for (size_t j = 0; j < c_; ++j)
                vector_[i][j] = vectors[i][j];
    }

//    void deallocate() {}

public:
    static Matrix zeros(size_t r, size_t c) { return Matrix(r, c); }

    static Matrix zeros(size_t n) { return Matrix(n); }

    static Matrix ones(size_t r, size_t c) { return Matrix(r, c, 1); }

    static Matrix ones(size_t n) { return Matrix(n, n, 1); }

    static Matrix identity(size_t n) { return diagonal(n, 1); }

    static Matrix id(size_t n) { return identity(n); }

    static Matrix diagonal(size_t n, T_ value) {
        Matrix tmp(n);
        for (size_t i = 0; i < n; ++i) tmp[i][i] = value;
        return tmp;
    }

    template<typename... Vectors>
    static Matrix from(const Vectors &...vectors) {
        Matrix tmp(sizeof...(vectors), 0);
        tmp.put_array(0, vectors...);
        return tmp;
    }

    template<size_t N, typename... Arrays>
    static Matrix from(const Arrays (&...arrays)[N]) {
        Matrix tmp(N, 0);
        tmp.put_array(0, arrays...);
        return tmp;
    }

    static Matrix from_row(const Vector<T_> &vector) { return Matrix::from(vector); }

    static Matrix from_col(const Vector<T_> &vector) {
        Matrix result(vector.size_, 1);
        for (size_t i = 0; i < vector.size_; ++i) result[i][0] = vector.arr_[i];
        return result;
    }

    static Matrix quad(const Matrix M) { return Matrix(M, M, M, M); }
};

template<typename T_>
Matrix<T_> operator*(T_ lhs, const Matrix<T_> &rhs) {
    Matrix<T_> tmp(rhs);
    tmp.operator*=(lhs);
    return tmp;
}

template<typename T_>
T_ det(const Matrix<T_> &A) { return A.det(); }

template<typename T_>
T_ tr(const Matrix<T_> &A) { return A.tr(); }

template<typename T_>
Matrix<T_> inv(const Matrix<T_> &A) { return A.inv(); }

template<typename T_>
Matrix<T_> RRE(const Matrix<T_> &A) { return A.RRE(); }

template<typename T_>
class MatrixLU {
private:
    using PM = pair<Matrix<T_>, Matrix<T_>>;
    Matrix<T_> l_;
    Matrix<T_> u_;

public:
    MatrixLU(const MatrixLU &other) : l_(other.l_), u_(other.u_) {}

    MatrixLU(MatrixLU &&other) noexcept: l_(move(other.l_)), u_(move(other.u_)) {}

    explicit MatrixLU(const PM &lu) : l_(lu.first), u_(lu.second) {}

    MatrixLU(const Matrix<T_> &l, const Matrix<T_> &u) : l_(l), u_(u) {}

    Matrix<T_> &l() { return l_; }

    constexpr const Matrix<T_> &l() const { return l_; }

    Matrix<T_> &u() { return u_; }

    constexpr const Matrix<T_> &u() const { return u_; }
};

using numeric_matrix = Matrix<double>;

#endif //VNET_MATRIX_H
