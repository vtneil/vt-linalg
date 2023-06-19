/**
 * @file dynamic_numeric_matrix.h
 * @author Vivatsathorn Thitasirivit
 * @date 31 May 2023
 * @brief Numeric matrix library
 */

#ifndef VNET_LINALG_DYNAMIC_NUMERIC_MATRIX_H
#define VNET_LINALG_DYNAMIC_NUMERIC_MATRIX_H

#include "standard_utility.h"
#include "pair.h"
#include "iterator.h"
#include "dynamic_numeric_vector.h"

namespace vt {
    template<typename T>
    class numeric_matrix_dynamic_lu_t;

    template<typename T>
    class numeric_matrix_dynamic_t {
    private:
        friend class numeric_vector_dynamic_t<T>;

    private:
        static constexpr size_t STRASSEN_DIMENSION = 1024;
        static constexpr size_t STRASSEN_THRESHOLD = STRASSEN_DIMENSION * STRASSEN_DIMENSION;
        size_t r_;
        size_t c_;
        numeric_vector_dynamic_t<numeric_vector_dynamic_t<T>> vector_;

    public:
        numeric_matrix_dynamic_t() : r_(0), c_(0), vector_(vt::move(numeric_vector_dynamic_t<numeric_vector_dynamic_t<T>>())) {}

        explicit numeric_matrix_dynamic_t(size_t n) : r_(n), c_(n) { allocate_zero(); }

        numeric_matrix_dynamic_t(size_t rows, size_t cols) : r_(rows), c_(cols) { allocate_zero(); }

        numeric_matrix_dynamic_t(size_t rows, size_t cols, T fill) : r_(rows), c_(cols) { allocate_fill(fill); }

        numeric_matrix_dynamic_t(const numeric_matrix_dynamic_t &other) : r_(other.r_), c_(other.c_) { allocate_from(other); }

        numeric_matrix_dynamic_t(numeric_matrix_dynamic_t &&other) noexcept: r_(other.r_), c_(other.c_) {
            steal(vt::move(other));
            other.r_ = 0;
            other.c_ = 0;
        }

        template<size_t R, size_t C>
        explicit numeric_matrix_dynamic_t(const T (&array)[R][C]) : r_(R), c_(C) { allocate_from(array); }

        template<size_t R>
        explicit numeric_matrix_dynamic_t(const numeric_vector_dynamic_t<T> (&vectors)[R]) : r_(R), c_(vectors[0].size_) { allocate_from(vectors); }

        numeric_matrix_dynamic_t(const numeric_matrix_dynamic_t &M11, const numeric_matrix_dynamic_t &M12, const numeric_matrix_dynamic_t &M21, const numeric_matrix_dynamic_t &M22) {
            if (M11.r_ == M12.r_ && M21.r_ == M22.r_ && M11.c_ == M21.c_ && M12.c_ == M22.c_) {
                r_ = M11.r_ + M21.r_;
                c_ = M11.c_ + M12.c_;
                allocate_zero();
                insert(0, 0, M11);
                insert(0, M11.c_, M12);
                insert(M11.r_, 0, M21);
                insert(M11.r_ / 2, M11.c_ / 2, M22);
            } else {
                r_ = vt::max(M11.r_ + M21.r_, M12.r_ + M22.r_);
                c_ = vt::max(M11.c_ + M12.c_, M21.c_ + M22.c_);
                allocate_zero();
            }
        }

        numeric_vector_dynamic_t<T> &operator[](size_t index) { return vector_[index]; }

        const numeric_vector_dynamic_t<T> &operator[](size_t index) const { return vector_[index]; }

        T &at(size_t r_index, size_t c_index) { return vector_[r_index][c_index]; };

        const T &at(size_t r_index, size_t c_index) const { return vector_[r_index][c_index]; };

        T &operator()(size_t r_index, size_t c_index) { return at(r_index, c_index); }

        const T &operator()(size_t r_index, size_t c_index) const { return at(r_index, c_index); }

        numeric_matrix_dynamic_t &operator=(const numeric_matrix_dynamic_t &other) {
            if (this != &other) {
                if (r_ == other.r_ && c_ == other.c_) {
                    insert(other);
                } else {
                    allocate_from(other);
                }
            }
            return *this;
        }

        numeric_matrix_dynamic_t &operator=(numeric_matrix_dynamic_t &&other) noexcept {
            if (this != &other) {
                steal(vt::move(other));
                other.r_ = 0;
                other.c_ = 0;
            }
            return *this;
        }

        template<size_t R, size_t C>
        numeric_matrix_dynamic_t &operator=(const T (&array)[R][C]) {
            if (r_ == R && c_ == C) {
                insert(array);
            } else {
                allocate_from(array);
            }
            return *this;
        }

        template<size_t R>
        numeric_matrix_dynamic_t &operator=(const numeric_vector_dynamic_t<T> (&vectors)[R]) {
            allocate_from(vectors);
            return *this;
        }

        numeric_matrix_dynamic_t &operator+=(const numeric_matrix_dynamic_t &other) { return iadd(*this, *this, other); }

        numeric_matrix_dynamic_t operator+(const numeric_matrix_dynamic_t &other) const {
            numeric_matrix_dynamic_t tmp(*this);
            tmp.operator+=(other);
            return tmp;
        }

        numeric_matrix_dynamic_t add(const numeric_matrix_dynamic_t &other) const { return operator+(other); }

        static numeric_matrix_dynamic_t &iadd(numeric_matrix_dynamic_t &C, const numeric_matrix_dynamic_t &A, const numeric_matrix_dynamic_t &B) {
            for (size_t i = 0; i < C.r_; ++i)
                for (size_t j = 0; j < C.c_; ++j)
                    C[i][j] = A[i][j] + B[i][j];
            return C;
        }

        numeric_matrix_dynamic_t &operator-=(const numeric_matrix_dynamic_t &other) { return isub(*this, *this, other); }

        numeric_matrix_dynamic_t operator-(const numeric_matrix_dynamic_t &other) const {
            numeric_matrix_dynamic_t tmp(*this);
            tmp.operator-=(other);
            return tmp;
        }

        numeric_matrix_dynamic_t sub(const numeric_matrix_dynamic_t &other) const { return operator-(other); }

        static numeric_matrix_dynamic_t &isub(numeric_matrix_dynamic_t &C, const numeric_matrix_dynamic_t &A, const numeric_matrix_dynamic_t &B) {
            for (size_t i = 0; i < C.r_; ++i)
                for (size_t j = 0; j < C.c_; ++j)
                    C[i][j] = A[i][j] - B[i][j];
            return C;
        }

        numeric_matrix_dynamic_t &operator*=(const numeric_matrix_dynamic_t &other) {
            steal(vt::move(operator*(other)));
            return *this;
        }

        numeric_matrix_dynamic_t operator*(const numeric_matrix_dynamic_t &other) const {
            numeric_matrix_dynamic_t C(r_, other.c_);
            return imatmul(C, *this, other);
        }

        numeric_matrix_dynamic_t &operator*=(T rhs) {
            for (auto &x: vector_) x *= rhs;
            return *this;
        }

        numeric_matrix_dynamic_t operator*(T rhs) const {
            numeric_matrix_dynamic_t tmp(*this);
            tmp.operator*=(rhs);
            return tmp;
        }

        numeric_vector_dynamic_t<T> operator*(const numeric_vector_dynamic_t<T> &other) const {
            numeric_vector_dynamic_t<T> tmp(r_);
            for (size_t i = 0; i < r_; ++i) tmp[i] = vector_[i].dot(other);
            return tmp;
        }

        numeric_matrix_dynamic_t matmul(const numeric_matrix_dynamic_t &other) const { return operator*(other); }

        static numeric_matrix_dynamic_t &imatmul(numeric_matrix_dynamic_t &C, const numeric_matrix_dynamic_t &A, const numeric_matrix_dynamic_t &B) {
            if (A.r() != 1 &&
                B.r() != 1 &&
                A.c() != 1 &&
                B.c() != 1 &&
                vt::max(A.r() * A.c(), B.r() * B.c()) >= STRASSEN_THRESHOLD)
                return mm_strassen(C, A, B);
            return mm_naive(C, A, B);
        }

        numeric_matrix_dynamic_t matmul_naive(const numeric_matrix_dynamic_t &other) const {
            numeric_matrix_dynamic_t C(r_, other.c_);
            return mm_naive(C, *this, other);
        }

        numeric_vector_dynamic_t<T> transform(const numeric_vector_dynamic_t<T> &other) const { return operator*(other); }

        numeric_matrix_dynamic_t &operator^=(size_t n) {
            steal(vt::move(operator^(n)));
            return *this;
        }

        numeric_matrix_dynamic_t operator^(size_t n) {
            numeric_matrix_dynamic_t product(*this);
            if (n == 0) return id(vt::min(r_, c_));
            if (n == 1) return product;
            product = vt::move(product.operator^(n / 2));
            if (n % 2 == 0) return product.operator*(product);
            else return product.operator*(product).operator*(*this);
        }

        numeric_matrix_dynamic_t matpow(size_t n) { return operator^(n); }

        numeric_matrix_dynamic_t matpow_naive(size_t n) {
            numeric_matrix_dynamic_t product = vt::move(id(r_));
            for (size_t i = 0; i < n; ++i) product *= *this;
            return product;
        }

        bool operator==(const numeric_matrix_dynamic_t &other) const {
            if (this == &other) return true;
            if (r_ != other.r_ || c_ != other.c_) return false;
            for (size_t i = 0; i < r_; ++i)
                if (vector_[i] != other.vector_[i])
                    return false;
            return true;
        }

        template<size_t R, size_t C>
        bool operator==(const T (&array)[R][C]) const {
            if (r_ != R || c_ != C) return false;
            for (size_t i = 0; i < r_; ++i)
                if (vector_[i] != array[i])
                    return false;
            return true;
        }

        bool operator!=(const numeric_matrix_dynamic_t &other) const { return !operator==(other); }

        template<size_t R, size_t C>
        bool operator!=(const T (&array)[R][C]) const { return !operator==(array); }

        bool equals(const numeric_matrix_dynamic_t &other) const { return operator==(other); }

        bool float_equals(const numeric_matrix_dynamic_t &other) const {
            if (this == &other) return true;
            if (r_ != other.r_ || c_ != other.c_) return false;
            for (size_t i = 0; i < r_; ++i)
                if (!vector_[i].float_equals(other.vector_[i]))
                    return false;
            return true;
        }

        template<size_t R, size_t C>
        bool equals(const T (&array)[R][C]) const { return operator==(array); }

        bool is_square() const { return r_ == c_; }

        bool can_multiply_with(const numeric_matrix_dynamic_t &other) const { return c_ == other.r_; }

        bool can_multiply_with(const numeric_vector_dynamic_t<T> &other) const { return c_ == other.size_; }

        numeric_vector_dynamic_t<T> row(size_t r_index) const { return numeric_vector_dynamic_t<T>(operator[](r_index)); }

        numeric_vector_dynamic_t<T> col(size_t c_index) const {
            numeric_vector_dynamic_t<T> tmp(r_);
            for (size_t i = 0; i < r_; ++i) tmp[i] = vector_[i][c_index];
            return tmp;
        }

        numeric_vector_dynamic_t<T> diag() {
            size_t min_dim = vt::min(r_, c_);
            numeric_vector_dynamic_t<T> tmp(min_dim);
            for (size_t i = 0; i < min_dim; ++i) tmp[i] = vector_[i][i];
            return tmp;
        }

        numeric_matrix_dynamic_t transpose() const {
            numeric_matrix_dynamic_t tmp(c_, r_);
            for (size_t i = 0; i < r_; ++i)
                for (size_t j = 0; j < c_; ++j)
                    tmp.vector_[j][i] = vector_[i][j];
            return tmp;
        }

        T det() const {
            numeric_matrix_dynamic_lu_t<T> lu = vt::move(LU());
            T det = 1;
            size_t r_swaps = 0;
            size_t n = vt::min(r_, c_);
            for (size_t i = 0; i < n; ++i) {
                det *= lu.u()[i][i];
                if (lu.u()[i][i] == 0) return 0;
                if (lu.l()[i][i] != 1) ++r_swaps;
            }
            return (r_swaps % 2 == 0) ? det : -det;
        }

        T tr() const {
            T acc = 0;
            size_t n = vt::min(r_, c_);
            for (size_t i = 0; i < n; ++i) acc += vector_[i][i];
            return acc;
        }

        numeric_matrix_dynamic_t inv() const {
            numeric_matrix_dynamic_lu_t<T> lu = vt::move(LU());
            return inv_ut(lu.u()) * inv_lt(lu.l());
        }

        constexpr numeric_matrix_dynamic_t inverse() const { return inv(); }

        numeric_matrix_dynamic_t RRE() const {
            numeric_matrix_dynamic_t m(*this);
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
                T val = m[r][lead];
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

        numeric_matrix_dynamic_lu_t<T> LU() const {
            size_t n = vt::min(r_, c_);
            numeric_matrix_dynamic_t lower(n, n);
            numeric_matrix_dynamic_t upper(n, n);
            for (size_t i = 0; i < n; ++i) {
                for (size_t k = 0; k < n; ++k) {
                    T sum_ = 0;
                    for (size_t j = 0; j < i; ++j) sum_ += lower[i][j] * upper[j][k];
                    upper[i][k] = vector_[i][k] - sum_;
                }
                for (size_t k = i; k < n; ++k) {
                    if (i == k) lower[i][i] = 1;
                    else {
                        T sum_ = 0;
                        for (size_t j = 0; j < i; ++j) sum_ += lower[k][j] * upper[j][i];
                        lower[k][i] = (vector_[k][i] - sum_) / upper[i][i];
                    }
                }
            }
            return {lower, upper};
        }

        void insert(const numeric_matrix_dynamic_t &M) { insert(0, 0, M); }

        template<size_t R, size_t C>
        void insert(const T (&array)[R][C]) { insert(0, 0, array); }

        void insert(size_t pos_row, size_t pos_col, const numeric_matrix_dynamic_t &M) {
            if (pos_row >= r_ || pos_col >= c_ || pos_row + M.r_ > r_ || pos_col + M.c_ > c_) return;
            for (size_t i = 0; i < M.r_; ++i)
                for (size_t j = 0; j < M.c_; ++j)
                    vector_[i + pos_row][j + pos_col] = M[i][j];
        }

        template<size_t R, size_t C>
        void insert(size_t pos_row, size_t pos_col, const T (&array)[R][C]) {
            if (pos_row >= r_ || pos_col >= c_ || pos_row + R > r_ || pos_col + C > c_) return;
            for (size_t i = 0; i < R; ++i)
                for (size_t j = 0; j < C; ++j)
                    vector_[i + pos_row][j + pos_col] = array[i][j];
        }

        numeric_matrix_dynamic_t slice(size_t r1, size_t c1, size_t r2, size_t c2) const {
            if (r1 > r2) vt::swap(r1, r2);
            if (c1 > c2) vt::swap(c1, c2);
            r2 = vt::min(r2, r_);
            c2 = vt::min(c2, c_);
            size_t r = r2 - r1;
            size_t c = c2 - c1;
            numeric_matrix_dynamic_t result(r, c);
            if (r == 0 || c == 0 || r1 >= r_ || c1 >= c_) return result;
            for (size_t i = 0; i < r; ++i)
                for (size_t j = 0; j < c; ++j)
                    result[i][j] = vector_[r1 + i][c1 + j];
            return result;
        }

        void swap(numeric_matrix_dynamic_t &other) {
            vt::swap(vector_, other.vector_);
            vt::swap(r_, other.r_);
            vt::swap(c_, other.c_);
        }

        iterator<numeric_vector_dynamic_t<T>> begin() { return vector_.begin(); }

        iterator<numeric_vector_dynamic_t<T>> begin() const { return vector_.begin(); }

        iterator<numeric_vector_dynamic_t<T>> end() { return vector_.end(); }

        iterator<numeric_vector_dynamic_t<T>> end() const { return vector_.end(); }

        constexpr size_t r() const { return r_; }

        constexpr size_t c() const { return c_; }

        constexpr size_t n() const { return r_ * c_; }

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
        static numeric_matrix_dynamic_t &inv_lt(numeric_matrix_dynamic_t &L) {
            size_t n = vt::min(L.r_, L.c_);
            numeric_matrix_dynamic_t X(n);
            for (size_t k = 0; k < n; ++k) {
                X[k][k] = 1 / L[k][k];
                for (size_t i = k + 1; i < n; ++i) {
                    T acc = 0;
                    for (size_t j = k; j < i; ++j)
                        acc -= L[i][j] * X[j][k];
                    X[i][k] = acc / L[i][i];
                }
            }
            L = vt::move(X);
            return L;
        }

        static numeric_matrix_dynamic_t inv_ut(const numeric_matrix_dynamic_t &U) {
            numeric_matrix_dynamic_t tmp = vt::move(U.transpose());
            return inv_lt(tmp).transpose();
        }

        static numeric_matrix_dynamic_t &mm_naive(numeric_matrix_dynamic_t &C, const numeric_matrix_dynamic_t &A, const numeric_matrix_dynamic_t &B) {
            for (size_t i = 0; i < A.r_; ++i) {
                numeric_vector_dynamic_t<T> &row_C = C[i];
                for (size_t k = 0; k < A.c_; ++k) {
                    const T &val_A = A[i][k];
                    for (size_t j = 0; j < B.c_; ++j) {
                        row_C[j] += val_A * B[k][j];
                    }
                }
            }
            return C;
        }

        static numeric_matrix_dynamic_t &mm_strassen(numeric_matrix_dynamic_t &C, const numeric_matrix_dynamic_t &A, const numeric_matrix_dynamic_t &B) {
            size_t rA = closest_2(A.r_);
            size_t cA = closest_2(A.c_);
            size_t rB = closest_2(B.r_);
            size_t cB = closest_2(B.c_);

            numeric_matrix_dynamic_t A11, A12, A21, A22, B11, B12, B21, B22;

            if (rA == A.r_ && cA == A.c_ && rB == B.r_ && cB == B.c_) {
                A11 = vt::move(A.slice(0, 0, rA / 2, cA / 2));
                A12 = vt::move(A.slice(0, cA / 2, rA / 2, cA));
                A21 = vt::move(A.slice(rA / 2, 0, rA, cA / 2));
                A22 = vt::move(A.slice(rA / 2, cA / 2, rA, cA));
                B11 = vt::move(B.slice(0, 0, rB / 2, cB / 2));
                B12 = vt::move(B.slice(0, cB / 2, rB / 2, cB));
                B21 = vt::move(B.slice(rB / 2, 0, rB, cB / 2));
                B22 = vt::move(B.slice(rB / 2, cB / 2, rB, cB));
            } else {
                numeric_matrix_dynamic_t Ap(rA, cA);
                Ap.insert(A);
                numeric_matrix_dynamic_t Bp(rB, cB);
                Bp.insert(B);
                A11 = vt::move(Ap.slice(0, 0, rA / 2, cA / 2));
                A12 = vt::move(Ap.slice(0, cA / 2, rA / 2, cA));
                A21 = vt::move(Ap.slice(rA / 2, 0, rA, cA / 2));
                A22 = vt::move(Ap.slice(rA / 2, cA / 2, rA, cA));
                B11 = vt::move(Bp.slice(0, 0, rB / 2, cB / 2));
                B12 = vt::move(Bp.slice(0, cB / 2, rB / 2, cB));
                B21 = vt::move(Bp.slice(rB / 2, 0, rB, cB / 2));
                B22 = vt::move(Bp.slice(rB / 2, cB / 2, rB, cB));
                C = vt::move(numeric_matrix_dynamic_t(rA, cB));
            }

            numeric_matrix_dynamic_t M1 = vt::move((A11 + A22) * (B11 + B22));
            numeric_matrix_dynamic_t M2 = vt::move((A21 + A22) * B11);
            numeric_matrix_dynamic_t M3 = vt::move(A11 * (B12 - B22));
            numeric_matrix_dynamic_t M4 = vt::move(A22 * (B21 - B11));
            numeric_matrix_dynamic_t M5 = vt::move((A11 + A12) * B22);
            numeric_matrix_dynamic_t M6 = vt::move((A21 - A11) * (B11 + B12));
            numeric_matrix_dynamic_t M7 = vt::move((A12 - A22) * (B21 + B22));
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

        void put_array(size_t index, const numeric_vector_dynamic_t<T> &vector) {
            vector_[index] = vector;
            c_ = vector.size();
        }

        template<size_t N>
        void put_array(size_t index, const T (&array)[N_]) {
            vector_[index] = array;
            c_ = N_;
        }

        template<typename... Vectors>
        void put_array(size_t index, const numeric_vector_dynamic_t<T> &vector, const Vectors &...vectors) {
            vector_[index] = vector;
            put_array(index + 1, vectors...);
        }

        template<size_t N, typename... Arrays>
        void put_array(size_t index, const T (&array)[N_], const Arrays (&...arrays)[N_]) {
            vector_[index] = array;
            put_array(index + 1, arrays...);
        }

        void allocate_zero() { vector_ = vt::move(numeric_vector_dynamic_t<numeric_vector_dynamic_t<T>>(r_, numeric_vector_dynamic_t<T>(c_))); }

        void allocate_fill(T fill) { vector_ = vt::move(numeric_vector_dynamic_t<numeric_vector_dynamic_t<T>>(r_, numeric_vector_dynamic_t<T>(c_, fill))); }

        void steal(numeric_matrix_dynamic_t &&other) {
            r_ = other.r_;
            c_ = other.c_;
            vector_ = vt::move(other.vector_);
        }

        void allocate_from(const numeric_matrix_dynamic_t &other) {
            r_ = other.r_;
            c_ = other.c_;
            vector_ = other.vector_;
        }

        template<size_t R, size_t C>
        void allocate_from(const T (&array)[R][C]) {
            r_ = R;
            c_ = C;
            allocate_zero();
            for (size_t i = 0; i < r_; ++i)
                for (size_t j = 0; j < c_; ++j)
                    vector_[i][j] = array[i][j];
        }

        template<size_t R>
        void allocate_from(const numeric_vector_dynamic_t<T> (&vectors)[R]) {
            r_ = R;
            c_ = vectors[0].size_;
            allocate_zero();
            for (size_t i = 0; i < r_; ++i)
                for (size_t j = 0; j < c_; ++j)
                    vector_[i][j] = vectors[i][j];
        }

        //    void deallocate() {}

    public:
        static numeric_matrix_dynamic_t zeros(size_t r, size_t c) { return numeric_matrix_dynamic_t(r, c); }

        static numeric_matrix_dynamic_t zeros(size_t n) { return numeric_matrix_dynamic_t(n); }

        static numeric_matrix_dynamic_t ones(size_t r, size_t c) { return numeric_matrix_dynamic_t(r, c, 1); }

        static numeric_matrix_dynamic_t ones(size_t n) { return numeric_matrix_dynamic_t(n, n, 1); }

        static numeric_matrix_dynamic_t identity(size_t n) { return diagonal(n, 1); }

        static numeric_matrix_dynamic_t id(size_t n) { return identity(n); }

        static numeric_matrix_dynamic_t diagonal(size_t n, T value) {
            numeric_matrix_dynamic_t tmp(n);
            for (size_t i = 0; i < n; ++i) tmp[i][i] = value;
            return tmp;
        }

        template<typename... Vectors>
        static numeric_matrix_dynamic_t from(const Vectors &...vectors) {
            numeric_matrix_dynamic_t tmp(sizeof...(vectors), 0);
            tmp.put_array(0, vectors...);
            return tmp;
        }

        template<size_t N, typename... Arrays>
        static numeric_matrix_dynamic_t from(const Arrays (&...arrays)[N_]) {
            numeric_matrix_dynamic_t tmp(N_, 0);
            tmp.put_array(0, arrays...);
            return tmp;
        }

        static numeric_matrix_dynamic_t from_row(const numeric_vector_dynamic_t<T> &vector) { return numeric_matrix_dynamic_t::from(vector); }

        static numeric_matrix_dynamic_t from_col(const numeric_vector_dynamic_t<T> &vector) {
            numeric_matrix_dynamic_t result(vector.size_, 1);
            for (size_t i = 0; i < vector.size_; ++i) result[i][0] = vector.arr_[i];
            return result;
        }

        static numeric_matrix_dynamic_t quad(const numeric_matrix_dynamic_t M) { return numeric_matrix_dynamic_t(M, M, M, M); }
    };

    template<typename T>
    numeric_matrix_dynamic_t<T> operator*(T lhs, const numeric_matrix_dynamic_t<T> &rhs) {
        numeric_matrix_dynamic_t < T > tmp(rhs);
        tmp.operator*=(lhs);
        return tmp;
    }

    template<typename T>
    T det(const numeric_matrix_dynamic_t<T> &A) { return A.det(); }

    template<typename T>
    T tr(const numeric_matrix_dynamic_t<T> &A) { return A.tr(); }

    template<typename T>
    numeric_matrix_dynamic_t<T> inv(const numeric_matrix_dynamic_t<T> &A) { return A.inv(); }

    template<typename T>
    numeric_matrix_dynamic_t<T> RRE(const numeric_matrix_dynamic_t<T> &A) { return A.RRE(); }

    template<typename T>
    class numeric_matrix_dynamic_lu_t {
    private:
        using PM = vt::pair<numeric_matrix_dynamic_t<T>, numeric_matrix_dynamic_t<T>>;
        numeric_matrix_dynamic_t<T> l_;
        numeric_matrix_dynamic_t<T> u_;

    public:
        numeric_matrix_dynamic_lu_t(const numeric_matrix_dynamic_lu_t &other) : l_(other.l_), u_(other.u_) {}

        numeric_matrix_dynamic_lu_t(numeric_matrix_dynamic_lu_t &&other) noexcept: l_(vt::move(other.l_)), u_(vt::move(other.u_)) {}

        explicit numeric_matrix_dynamic_lu_t(const PM &lu) : l_(lu.first), u_(lu.second) {}

        numeric_matrix_dynamic_lu_t(const numeric_matrix_dynamic_t<T> &l, const numeric_matrix_dynamic_t<T> &u) : l_(l), u_(u) {}

        numeric_matrix_dynamic_t<T> &l() { return l_; }

        constexpr const numeric_matrix_dynamic_t<T> &l() const { return l_; }

        numeric_matrix_dynamic_t<T> &u() { return u_; }

        constexpr const numeric_matrix_dynamic_t<T> &u() const { return u_; }
    };

    using numeric_matrix = numeric_matrix_dynamic_t<double>;
}

#endif //VNET_LINALG_DYNAMIC_NUMERIC_MATRIX_H
