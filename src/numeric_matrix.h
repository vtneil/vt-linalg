/**
 * @file numeric_matrix.h
 * @author Vivatsathorn Thitasirivit
 * @date 15 June 2023
 * @brief Numeric matrix library
 */

#ifndef VNET_LINALG_NUMERIC_MATRIX_H
#define VNET_LINALG_NUMERIC_MATRIX_H

#include "standard_utility.h"
#include "pair.h"
#include "iterator.h"
#include "numeric_vector.h"

namespace vt {
    template<typename T, size_t OSize>
    class numeric_matrix_static_lu_t;

    template<typename T, size_t Row, size_t Col>
    class numeric_matrix_static_t {
    public:
        static_assert(Row > 0, "Row must be greater than 0.");
        static_assert(Col > 0, "Column must be greater than 0.");

    private:
        template<typename U, size_t V>
        friend
        class numeric_vector_static_t;

    private:
        static constexpr size_t Order = (Row < Col) ? Row : Col;
        numeric_vector_static_t<numeric_vector_static_t<T, Col>, Row> vector_ = {};

    public:
        constexpr numeric_matrix_static_t() = default;

        explicit numeric_matrix_static_t(T fill) { allocate_fill(fill); }

        numeric_matrix_static_t(const numeric_matrix_static_t &other) { allocate_from(other); }

        numeric_matrix_static_t(numeric_matrix_static_t &&other) noexcept { steal(vt::move(other)); }

        explicit numeric_matrix_static_t(const T (&array)[Row][Col]) { allocate_from(array); }

        explicit numeric_matrix_static_t(const numeric_vector_static_t<T, Col> (&vectors)[Row]) {
            allocate_from(vectors);
        }

        template<size_t R1, size_t R2, size_t C1, size_t C2>
        numeric_matrix_static_t(const numeric_matrix_static_t<T, R1, C1> &M11,
                                const numeric_matrix_static_t<T, R1, C2> &M12,
                                const numeric_matrix_static_t<T, R2, C1> &M21,
                                const numeric_matrix_static_t<T, R2, C2> &M22) {
            static_assert(R1 + R2 == Row, "Row must be the sum of R1 and R2.");
            static_assert(C1 + C2 == Col, "Column must be the sum of C1 and C2.");
            insert<0, 0>(M11);
            insert<0, C1>(M12);
            insert<R1, 0>(M21);
            insert<R1, C1>(M22);
        }

        numeric_vector_static_t<T, Col> &operator[](size_t index) { return vector_[index]; }

        constexpr const numeric_vector_static_t<T, Col> &operator[](size_t index) const { return vector_[index]; }

        T &at(size_t r_index, size_t c_index) { return vector_[r_index][c_index]; };

        constexpr const T &at(size_t r_index, size_t c_index) const { return vector_[r_index][c_index]; };

        T &operator()(size_t r_index, size_t c_index) { return at(r_index, c_index); }

        constexpr const T &operator()(size_t r_index, size_t c_index) const { return at(r_index, c_index); }

        numeric_matrix_static_t &operator=(const numeric_matrix_static_t &other) {
            if (this != &other) allocate_from(other);
            return *this;
        }

        numeric_matrix_static_t &operator=(numeric_matrix_static_t &&other) noexcept {
            if (this != &other) steal(vt::move(other));
            return *this;
        }

        numeric_matrix_static_t &operator=(const T (&array)[Row][Col]) {
            allocate_from(array);
            return *this;
        }

        numeric_matrix_static_t &operator=(const numeric_vector_static_t<T, Col> (&vectors)[Row]) {
            allocate_from(vectors);
            return *this;
        }

        numeric_matrix_static_t &operator+=(const numeric_matrix_static_t &other) { return iadd(*this, *this, other); }

        numeric_matrix_static_t operator+(const numeric_matrix_static_t &other) const {
            numeric_matrix_static_t tmp(*this);
            tmp.operator+=(other);
            return tmp;
        }

        constexpr numeric_matrix_static_t add(const numeric_matrix_static_t &other) const { return operator+(other); }

        static numeric_matrix_static_t &
        iadd(numeric_matrix_static_t &C, const numeric_matrix_static_t &A, const numeric_matrix_static_t &B) {
            for (size_t i = 0; i < Row; ++i)
                for (size_t j = 0; j < Col; ++j)
                    C[i][j] = A[i][j] + B[i][j];
            return C;
        }

        numeric_matrix_static_t &operator-=(const numeric_matrix_static_t &other) { return isub(*this, *this, other); }

        numeric_matrix_static_t operator-(const numeric_matrix_static_t &other) const {
            numeric_matrix_static_t tmp(*this);
            tmp.operator-=(other);
            return tmp;
        }

        constexpr numeric_matrix_static_t sub(const numeric_matrix_static_t &other) const { return operator-(other); }

        static numeric_matrix_static_t &
        isub(numeric_matrix_static_t &C, const numeric_matrix_static_t &A, const numeric_matrix_static_t &B) {
            for (size_t i = 0; i < Row; ++i)
                for (size_t j = 0; j < Col; ++j)
                    C[i][j] = A[i][j] - B[i][j];
            return C;
        }

        numeric_matrix_static_t<T, Order, Order> &operator*=(const numeric_matrix_static_t<T, Order, Order> &other) {
            steal(vt::move(operator*(other)));
            return *this;
        }

        template<size_t ORow, size_t OCol>
        numeric_matrix_static_t<T, Row, OCol> operator*(const numeric_matrix_static_t<T, ORow, OCol> &other) const {
            numeric_matrix_static_t<T, Row, OCol> C;
            return imatmul(C, *this, other);
        }

        numeric_matrix_static_t &operator*=(T rhs) {
            for (auto &x: vector_) x *= rhs;
            return *this;
        }

        numeric_matrix_static_t operator*(T rhs) const {
            numeric_matrix_static_t tmp(*this);
            tmp.operator*=(rhs);
            return tmp;
        }

        numeric_vector_static_t<T, Row> operator*(const numeric_vector_static_t<T, Col> &other) const {
            numeric_vector_static_t<T, Row> tmp;
            for (size_t i = 0; i < Row; ++i) tmp[i] = vector_[i].dot(other);
            return tmp;
        }

        template<size_t ORow, size_t X, size_t OCol>
        static constexpr numeric_matrix_static_t<T, ORow, OCol> imatmul(numeric_matrix_static_t<T, ORow, OCol> &C,
                                                                        const numeric_matrix_static_t<T, ORow, X> &A,
                                                                        const numeric_matrix_static_t<T, X, OCol> &B) {
            return mm_naive(C, A, B);
        }

        template<size_t ORow, size_t OCol>
        numeric_matrix_static_t<T, Row, OCol> matmul_naive(const numeric_matrix_static_t<T, ORow, OCol> &other) const {
            numeric_matrix_static_t<T, Row, OCol> C;
            return mm_naive(C, *this, other);
        }

        constexpr numeric_vector_static_t<T, Row> transform(const numeric_vector_static_t<T, Col> &other) const {
            return operator*(other);
        }

        numeric_matrix_static_t<T, Order, Order> &operator^=(size_t n) {
            static_assert(static_is_a_square_matrix(), "Non-square matrix can\'t use power operator.");
            steal(vt::move(operator^(n)));
            return *this;
        }

        numeric_matrix_static_t<T, Order, Order> operator^(size_t n) {
            static_assert(static_is_a_square_matrix(), "Non-square matrix can\'t use power operator.");
            numeric_matrix_static_t<T, Order, Order> product(*this);
            if (n == 0) return id();
            if (n == 1) return product;
            product = vt::move(product.operator^(n / 2));
            if (n % 2 == 0) return product.operator*(product);
            else return product.operator*(product).operator*(*this);
        }

        numeric_matrix_static_t<T, Order, Order> matpow(size_t n) {
            static_assert(static_is_a_square_matrix(), "Non-square matrix can\'t use power operator.");
            return operator^(n);
        }

        numeric_matrix_static_t<T, Order, Order> matpow_naive(size_t n) {
            static_assert(static_is_a_square_matrix(), "Non-square matrix can\'t use power operator.");
            numeric_matrix_static_t<T, Order, Order> product = vt::move(id());
            for (size_t i = 0; i < n; ++i) product *= *this;
            return product;
        }

        template<size_t pos_row = 0, size_t pos_col = 0, size_t ORow, size_t OCol>
        void insert(const numeric_matrix_static_t<T, ORow, OCol> &M) {
            static_assert(pos_row < Row, "Insertion failed! Start row is out of range.");
            static_assert(pos_col < Col, "Insertion failed! Start column is out of range.");
            static_assert(pos_row + ORow <= Row, "Insertion failed! Stop row is out of range.");
            static_assert(pos_col + OCol <= Col, "Insertion failed! Stop column is out of range.");
            for (size_t i = 0; i < ORow; ++i)
                for (size_t j = 0; j < OCol; ++j)
                    vector_[pos_row + i][pos_col + j] = M[i][j];
        }

        template<size_t pos_row = 0, size_t pos_col = 0, size_t ORow, size_t OCol>
        void insert(const T (&array)[ORow][OCol]) {
            static_assert(pos_row < Row, "Insertion failed! Start row is out of range.");
            static_assert(pos_col < Col, "Insertion failed! Start column is out of range.");
            static_assert(pos_row + ORow <= Row, "Insertion failed! Stop row is out of range.");
            static_assert(pos_col + OCol <= Col, "Insertion failed! Stop column is out of range.");
            for (size_t i = 0; i < ORow; ++i)
                for (size_t j = 0; j < OCol; ++j)
                    vector_[pos_row + i][pos_col + j] = array[i][j];
        }

        template<size_t r1, size_t c1, size_t r2, size_t c2>
        numeric_matrix_static_t<T, r2 - r1, c2 - c1> slice() const {
            static_assert(r1 < r2, "Start row must be less than stop row.");
            static_assert(c1 < c2, "Start column must be less than stop column.");
            static_assert(r2 <= Row, "Row is out of range.");
            static_assert(c2 <= Col, "Column is out of range.");
            numeric_matrix_static_t<T, r2 - r1, c2 - c1> result;
            for (size_t i = 0; i < r2 - r1; ++i)
                for (size_t j = 0; j < c2 - c1; ++j)
                    result[i][j] = vector_[r1 + i][c1 + j];
            return result;
        }

        constexpr numeric_vector_static_t<T, Col> row(size_t r_index) const {
            return numeric_vector_static_t<T, Col>(operator[](r_index));
        }

        numeric_vector_static_t<T, Row> col(size_t c_index) const {
            numeric_vector_static_t<T, Row> result;
            for (size_t i = 0; i < Row; ++i) result[i] = vector_[i][c_index];
            return result;
        }

        numeric_vector_static_t<T, Order> diag() const {
            numeric_vector_static_t<T, Order> result;
            for (size_t i = 0; i < Order; ++i) result[i] = vector_[i][i];
            return result;
        }

        numeric_matrix_static_t<T, Col, Row> transpose() const {
            numeric_matrix_static_t<T, Col, Row> result;
            for (size_t i = 0; i < Row; ++i)
                for (size_t j = 0; j < Col; ++j)
                    result[j][i] = vector_[i][j];
            return result;
        }

        T det() const {
            static_assert(static_is_a_square_matrix(), "Can only find determinant of a square matrix.");
            numeric_matrix_static_lu_t<T, Order> lu = vt::move(LU());
            T det = 1;
            size_t r_swaps = 0;
            for (size_t i = 0; i < Order; ++i) {
                det *= lu.u()[i][i];
                if (lu.u()[i][i] == 0) return 0;
                if (lu.l()[i][i] != 1) r_swaps++;
            }
            return (r_swaps % 2 == 0) ? det : -det;
        }

        T tr() const {
            static_assert(static_is_a_square_matrix(), "Can only find trace of a square matrix.");
            T acc = 0;
            for (size_t i = 0; i < Order; ++i) acc += vector_[i][i];
            return acc;
        }

        numeric_matrix_static_t inv() const {
            static_assert(static_is_a_square_matrix(), "Can only find inverse of a square matrix.");
            numeric_matrix_static_lu_t<T, Order> lu = vt::move(LU());
            return inv_ut(lu.u()) * inv_lt(lu.l());
        }

        constexpr numeric_matrix_static_t inverse() const { return inv(); }

        numeric_matrix_static_lu_t<T, Order> LU() const {
            static_assert(static_is_a_square_matrix(), "Can only find LU decomposition of a square matrix.");
            numeric_matrix_static_t<T, Order, Order> lower;
            numeric_matrix_static_t<T, Order, Order> upper;
            for (size_t i = 0; i < Order; ++i) {
                for (size_t k = 0; k < Order; ++k) {
                    T sum_ = 0;
                    for (size_t j = 0; j < i; ++j) sum_ += lower[i][j] * upper[j][k];
                    upper[i][k] = vector_[i][k] - sum_;
                }
                for (size_t k = i; k < Order; ++k) {
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

        numeric_matrix_static_t RRE() const {
            numeric_matrix_static_t m(*this);
            size_t lead = 0;
            for (size_t r = 0; r < Row; ++r) {
                if (lead >= Col) return m;
                size_t i;
                for (i = r; m[i][lead] == 0;) {
                    ++i;
                    if (i == Row) {
                        i = r;
                        ++lead;
                        if (lead == Col) return m;
                    }
                }
                m[i].swap(m[r]);
                T val = m[r][lead];
                for (size_t j = 0; j < Col; ++j) m[r][j] /= val;
                for (i = 0; i < Row; ++i) {
                    if (i != r) {
                        val = m[i][lead];
                        for (size_t j = 0; j < Col; ++j)
                            m[i][j] -= val * m[r][j];
                    }
                }
                ++lead;
            }
            m.fix_zero();
            return m;
        }

        template<size_t ORow, size_t OCol>
        bool operator==(const numeric_matrix_static_t<T, ORow, OCol> &other) const {
            if (this == &other) return true;
            if (Row != ORow || Col != OCol) return false;
            for (size_t i = 0; i < Row; ++i) if (vector_[i] != other.vector_[i]) return false;
            return true;
        }

        template<size_t ORow, size_t OCol>
        bool operator==(const T (&array)[ORow][OCol]) const {
            if (Row != ORow || Col != OCol) return false;
            for (size_t i = 0; i < Row; ++i) if (vector_[i] != array[i]) return false;
            return true;
        }

        template<size_t ORow, size_t OCol>
        constexpr bool operator!=(const numeric_matrix_static_t<T, ORow, OCol> &other) const {
            return !operator==(other);
        }

        template<size_t ORow, size_t OCol>
        constexpr bool operator!=(const T (&array)[ORow][OCol]) const { return !operator==(array); }

        template<size_t ORow, size_t OCol>
        constexpr bool equals(const numeric_matrix_static_t<T, ORow, OCol> &other) const { return operator==(other); }

        template<size_t ORow, size_t OCol>
        bool float_equals(const numeric_matrix_static_t<T, ORow, OCol> &other) const {
            if (this == &other) return true;
            if (Row != ORow || Col != OCol) return false;
            for (size_t i = 0; i < Row; ++i) if (!vector_[i].float_equals(other.vector_[i])) return false;
            return true;
        }

        template<size_t ORow, size_t OCol>
        constexpr bool equals(const T (&array)[ORow][OCol]) const { return operator==(array); }

        iterator<numeric_vector_static_t<T, Col>> begin() { return vector_.begin(); }

        constexpr iterator<numeric_vector_static_t<T, Col>> begin() const { return vector_.begin(); }

        iterator<numeric_vector_static_t<T, Col>> end() { return vector_.end(); }

        constexpr iterator<numeric_vector_static_t<T, Col>> end() const { return vector_.end(); }

        constexpr size_t r() const { return Row; }

        constexpr size_t c() const { return Col; }

        constexpr size_t n() const { return Row * Col; }

        constexpr size_t order() const { return Order; }

        constexpr bool is_square() const { return Row == Col; }

        void swap(numeric_matrix_static_t &other) {
            for (size_t i = 0; i < Row; ++i) vector_[i].swap(other.vector_[i]);
        }

        constexpr numeric_matrix_static_t copy() const { return numeric_matrix_static_t(*this); }

    private:
        void fix_zero() {
            for (size_t i = 0; i < Row; ++i)
                for (size_t j = 0; j < Col; ++j)
                    if (vector_[i][j] == 0.0)
                        vector_[i][j] = 0.0;
        }

        void put_array(size_t index, const numeric_vector_static_t<T, Col> &vector) { vector_[index] = vector; }

        void put_array(size_t index, const T (&array)[Col]) { vector_[index] = array; }

        template<typename... Vectors>
        void put_array(size_t index, const numeric_vector_static_t<T, Col> &vector, const Vectors &...vectors) {
            vector_[index] = vector;
            if (index < Row - 1) put_array(index + 1, vectors...);
        }

        template<typename... Arrays>
        void put_array(size_t index, const T (&array)[Col], const Arrays (&...arrays)[Col]) {
            vector_[index] = array;
            if (index < Row - 1) put_array(index + 1, arrays...);
        }

        void allocate_zero() { vector_ = vt::move(numeric_vector_static_t<numeric_vector_static_t<T, Col>, Row>()); }

        void allocate_fill(T fill) {
            vector_ = vt::move(
                    numeric_vector_static_t<numeric_vector_static_t<T, Col>, Row>(
                            numeric_vector_static_t<T, Col>(fill)));
        }

        void allocate_from(const numeric_matrix_static_t &other) { vector_ = other.vector_; }

        void allocate_from(const T (&array)[Row][Col]) {
            for (size_t i = 0; i < Row; ++i)
                for (size_t j = 0; j < Col; ++j)
                    vector_[i][j] = array[i][j];
        }

        void allocate_from(const numeric_vector_static_t<T, Col> (&vectors)[Row]) {
            for (size_t i = 0; i < Row; ++i)
                for (size_t j = 0; j < Col; ++j)
                    vector_[i][j] = vectors[i][j];
        }

        void steal(numeric_matrix_static_t &&other) { vector_ = vt::move(other.vector_); }

        static constexpr bool static_is_a_square_matrix() { return Row == Col; }

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

        template<size_t ORow, size_t X, size_t OCol>
        static numeric_matrix_static_t<T, ORow, OCol> &mm_naive(numeric_matrix_static_t<T, ORow, OCol> &C,
                                                                const numeric_matrix_static_t<T, ORow, X> &A,
                                                                const numeric_matrix_static_t<T, X, OCol> &B) {
            for (size_t i = 0; i < ORow; ++i) {
                numeric_vector_static_t<T, OCol> &row_C = C[i];
                for (size_t k = 0; k < X; ++k) {
                    const T &val_A = A[i][k];
                    for (size_t j = 0; j < OCol; ++j) {
                        row_C[j] += val_A * B[k][j];
                    }
                }
            }
            return C;
        }

        template<size_t ORow, size_t X, size_t OCol>
        static constexpr numeric_matrix_static_t<T, ORow, OCol> &mm_strassen(numeric_matrix_static_t<T, ORow, OCol> &C,
                                                                             const numeric_matrix_static_t<T, ORow, X> &A,
                                                                             const numeric_matrix_static_t<T, X, OCol> &B) {
            // to be implemented later.
            // not suitable for embedded system.
            return mm_naive(C, A, B);
        }

        template<size_t OSize>
        static numeric_matrix_static_t<T, OSize, OSize> inv_lt(numeric_matrix_static_t<T, OSize, OSize> &L) {
            static_assert(static_is_a_square_matrix(), "Can only inverse a square matrix.");
            numeric_matrix_static_t<T, OSize, OSize> X;
            for (size_t k = 0; k < OSize; ++k) {
                X[k][k] = 1 / L[k][k];
                for (size_t i = k + 1; i < OSize; ++i) {
                    T acc = 0;
                    for (size_t j = k; j < i; ++j)
                        acc -= L[i][j] * X[j][k];
                    X[i][k] = acc / L[i][i];
                }
            }
            L = vt::move(X);
            return L;
        }

        template<size_t OSize>
        static numeric_matrix_static_t<T, OSize, OSize> inv_ut(numeric_matrix_static_t<T, OSize, OSize> &U) {
            static_assert(static_is_a_square_matrix(), "Can only inverse a square matrix.");
            numeric_matrix_static_t<T, OSize, OSize> tmp = vt::move(U.transpose());
            return inv_lt(tmp).transpose();
        }

    public:
        static constexpr numeric_matrix_static_t zeros() { return numeric_matrix_static_t(); }

        static constexpr numeric_matrix_static_t ones() { return numeric_matrix_static_t(1); }

        static constexpr numeric_matrix_static_t identity() {
            static_assert(static_is_a_square_matrix(), "Identity matrix must be a square matrix.");
            return diagonal(1);
        }

        static constexpr numeric_matrix_static_t id() { return identity(); }

        static numeric_matrix_static_t diagonal(T value) {
            numeric_matrix_static_t tmp;
            for (size_t i = 0; i < Order; ++i) tmp[i][i] = value;
            return tmp;
        }

        template<size_t ORow, size_t OCol>
        static numeric_matrix_static_t constexpr quad(const numeric_matrix_static_t<T, ORow, OCol> M) {
            return numeric_matrix_static_t(M, M, M, M);
        }
    };

    template<typename T, size_t Row, size_t Col>
    numeric_matrix_static_t<T, Row, Col> operator*(T lhs, const numeric_matrix_static_t<T, Row, Col> &rhs) {
        numeric_matrix_static_t<T, Row, Col> tmp(rhs);
        tmp.operator*=(lhs);
        return tmp;
    }

    template<typename T, size_t Row, size_t Col>
    constexpr T det(const numeric_matrix_static_t<T, Row, Col> &A) { return A.det(); }

    template<typename T, size_t Row, size_t Col>
    constexpr T tr(const numeric_matrix_static_t<T, Row, Col> &A) { return A.tr(); }

    template<typename T, size_t Row, size_t Col>
    constexpr numeric_matrix_static_t<T, Row, Col> inv(const numeric_matrix_static_t<T, Row, Col> &A) {
        return A.inv();
    }

    template<typename T, size_t Row, size_t Col>
    constexpr numeric_matrix_static_t<T, Row, Col> RRE(const numeric_matrix_static_t<T, Row, Col> &A) {
        return A.RRE();
    }

    template<typename T, size_t OSize>
    class numeric_matrix_static_lu_t {
    private:
        using Matrix_t = numeric_matrix_static_t<T, OSize, OSize>;
        using Pair_t = vt::pair<numeric_matrix_static_t<T, OSize, OSize>, numeric_matrix_static_t<T, OSize, OSize>>;
        numeric_matrix_static_t<T, OSize, OSize> l_;
        numeric_matrix_static_t<T, OSize, OSize> u_;

    public:
        numeric_matrix_static_lu_t(const numeric_matrix_static_lu_t &other) : l_(other.l_), u_(other.u_) {}

        numeric_matrix_static_lu_t(numeric_matrix_static_lu_t &&other) noexcept: l_(vt::move(other.l_)),
                                                                                 u_(vt::move(other.u_)) {}

        explicit numeric_matrix_static_lu_t(const Pair_t &lu) : l_(lu.first), u_(lu.second) {}

        numeric_matrix_static_lu_t(const Matrix_t &l, const Matrix_t &u) : l_(l), u_(u) {}

        Matrix_t &l() { return l_; }

        constexpr const Matrix_t &l() const { return l_; }

        Matrix_t &u() { return u_; }

        constexpr const Matrix_t &u() const { return u_; }
    };

    template<size_t Row, size_t Col = Row>
    using numeric_matrix = numeric_matrix_static_t<real_t, Row, Col>;

    template<size_t Row, size_t Col = Row>
    constexpr numeric_matrix<Row, Col> make_numeric_matrix(const real_t (&array)[Row][Col]) {
        return numeric_matrix<Row, Col>(array);
    }
}

#endif //VNET_LINALG_NUMERIC_MATRIX_H
