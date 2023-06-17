/**
 * @file numeric_matrix.h
 * @author Vivatsathorn Thitasirivit
 * @date 15 June 2023
 * @brief Numeric matrix library
 */

#ifndef VNET_LINALG_STATIC_NUMERIC_MATRIX_H
#define VNET_LINALG_STATIC_NUMERIC_MATRIX_H

#include "utils.h"
#include "static_numeric_vector.h"

template<typename T, size_t OSize>
class MatrixStaticLU;

template<typename T, size_t Row, size_t Col>
class MatrixStatic {
public:
    static_assert(Row > 0, "Row must be greater than 0.");
    static_assert(Col > 0, "Column must be greater than 0.");

private:
    template<typename U, size_t V>
    friend
    class VectorStatic;

private:
    static constexpr size_t STRASSEN_DIMENSION = 1024;
    static constexpr size_t STRASSEN_THRESHOLD = STRASSEN_DIMENSION * STRASSEN_DIMENSION;
    static constexpr size_t Order = (Row < Col) ? Row : Col;
    VectorStatic<VectorStatic<T, Col>, Row> vector_ = {};

public:
    MatrixStatic() = default;

    explicit MatrixStatic(T fill) { allocate_fill(fill); }

    MatrixStatic(const MatrixStatic &other) { allocate_from(other); }

    MatrixStatic(MatrixStatic &&other) noexcept { steal(vt::move(other)); }

    explicit MatrixStatic(const T (&array)[Row][Col]) { allocate_from(array); }

    explicit MatrixStatic(const VectorStatic<T, Col> (&vectors)[Row]) { allocate_from(vectors); }

    template<size_t R1, size_t R2, size_t C1, size_t C2>
    MatrixStatic(const MatrixStatic<T, R1, C1> &M11,
                 const MatrixStatic<T, R1, C2> &M12,
                 const MatrixStatic<T, R2, C1> &M21,
                 const MatrixStatic<T, R2, C2> &M22) {
        static_assert(R1 + R2 == Row, "Row must be the sum of R1 and R2.");
        static_assert(C1 + C2 == Col, "Column must be the sum of C1 and C2.");
        insert(0, 0, M11);
        insert(0, C1, M12);
        insert(R1, 0, M21);
        insert(R1, C1, M22);
    }

    VectorStatic<T, Col> &operator[](size_t index) { return vector_[index]; }

    const VectorStatic<T, Col> &operator[](size_t index) const { return vector_[index]; }

    T &at(size_t r_index, size_t c_index) { return vector_[r_index][c_index]; };

    const T &at(size_t r_index, size_t c_index) const { return vector_[r_index][c_index]; };

    T &operator()(size_t r_index, size_t c_index) { return at(r_index, c_index); }

    const T &operator()(size_t r_index, size_t c_index) const { return at(r_index, c_index); }

    MatrixStatic &operator=(const MatrixStatic &other) {
        if (this != &other) allocate_from(other);
        return *this;
    }

    MatrixStatic &operator=(MatrixStatic &&other) noexcept {
        if (this != &other) steal(vt::move(other));
        return *this;
    }

    MatrixStatic &operator=(const T (&array)[Row][Col]) {
        allocate_from(array);
        return *this;
    }

    MatrixStatic &operator=(const VectorStatic<T, Col> (&vectors)[Row]) {
        allocate_from(vectors);
        return *this;
    }

    MatrixStatic &operator+=(const MatrixStatic &other) { return iadd(*this, *this, other); }

    MatrixStatic operator+(const MatrixStatic &other) const {
        MatrixStatic tmp(*this);
        tmp.operator+=(other);
        return tmp;
    }

    MatrixStatic add(const MatrixStatic &other) const { return operator+(other); }

    static MatrixStatic &iadd(MatrixStatic &C, const MatrixStatic &A, const MatrixStatic &B) {
        for (size_t i = 0; i < Row; ++i)
            for (size_t j = 0; j < Col; ++j)
                C[i][j] = A[i][j] + B[i][j];
        return C;
    }

    MatrixStatic &operator-=(const MatrixStatic &other) { return isub(*this, *this, other); }

    MatrixStatic operator-(const MatrixStatic &other) const {
        MatrixStatic tmp(*this);
        tmp.operator-=(other);
        return tmp;
    }

    MatrixStatic sub(const MatrixStatic &other) const { return operator-(other); }

    static MatrixStatic &isub(MatrixStatic &C, const MatrixStatic &A, const MatrixStatic &B) {
        for (size_t i = 0; i < Row; ++i)
            for (size_t j = 0; j < Col; ++j)
                C[i][j] = A[i][j] - B[i][j];
        return C;
    }

    MatrixStatic<T, Order, Order> &operator*=(const MatrixStatic<T, Order, Order> &other) {
        steal(vt::move(operator*(other)));
        return *this;
    }

    template<size_t ORow, size_t OCol>
    MatrixStatic<T, Row, OCol> operator*(const MatrixStatic<T, ORow, OCol> &other) const {
        MatrixStatic<T, Row, OCol> C;
        return imatmul(C, *this, other);
    }

    MatrixStatic &operator*=(T rhs) {
        for (auto &x: vector_) x *= rhs;
        return *this;
    }

    MatrixStatic operator*(T rhs) const {
        MatrixStatic tmp(*this);
        tmp.operator*=(rhs);
        return tmp;
    }

    VectorStatic<T, Row> operator*(const VectorStatic<T, Col> &other) const {
        VectorStatic<T, Row> tmp;
        for (size_t i = 0; i < Row; ++i) tmp[i] = vector_[i].dot(other);
        return tmp;
    }

    template<size_t ORow, size_t X, size_t OCol>
    static MatrixStatic<T, ORow, OCol> imatmul(MatrixStatic<T, ORow, OCol> &C,
                                               const MatrixStatic<T, ORow, X> &A,
                                               const MatrixStatic<T, X, OCol> &B) {
        return mm_naive(C, A, B);
    }

    template<size_t ORow, size_t OCol>
    MatrixStatic<T, Row, OCol> matmul_naive(const MatrixStatic<T, ORow, OCol> &other) const {
        MatrixStatic<T, Row, OCol> C;
        return mm_naive(C, *this, other);
    }

    VectorStatic<T, Row> transform(const VectorStatic<T, Col> &other) const { return operator*(other); }

    MatrixStatic<T, Order, Order> &operator^=(size_t n) {
        static_assert(is_square_(), "Non-square matrix can\'t use power operator.");
        steal(vt::move(operator^(n)));
        return *this;
    }

    MatrixStatic<T, Order, Order> operator^(size_t n) {
        static_assert(is_square_(), "Non-square matrix can\'t use power operator.");
        MatrixStatic<T, Order, Order> product(*this);
        if (n == 0) return id();
        if (n == 1) return product;
        product = vt::move(product.operator^(n / 2));
        if (n % 2 == 0) return product.operator*(product);
        else return product.operator*(product).operator*(*this);
    }

    MatrixStatic<T, Order, Order> matpow(size_t n) {
        static_assert(is_square_(), "Non-square matrix can\'t use power operator.");
        return operator^(n);
    }

    MatrixStatic<T, Order, Order> matpow_naive(size_t n) {
        static_assert(is_square_(), "Non-square matrix can\'t use power operator.");
        MatrixStatic<T, Order, Order> product = vt::move(id());
        for (size_t i = 0; i < n; ++i) product *= *this;
        return product;
    }

    template<size_t pos_row = 0, size_t pos_col = 0, size_t ORow, size_t OCol>
    void insert(const MatrixStatic<T, ORow, OCol> &M) {
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
    MatrixStatic<T, r2 - r1, c2 - c1> slice() const {
        static_assert(r1 < r2, "Start row must be less than stop row.");
        static_assert(c1 < c2, "Start column must be less than stop column.");
        static_assert(r2 <= Row, "Row is out of range.");
        static_assert(c2 <= Col, "Column is out of range.");
        MatrixStatic<T, r2 - r1, c2 - c1> result;
        for (size_t i = 0; i < r2 - r1; ++i)
            for (size_t j = 0; j < c2 - c1; ++j)
                result[i][j] = vector_[r1 + i][c1 + j];
        return result;
    }

    VectorStatic<T, Col> row(size_t r_index) const { return VectorStatic<T, Col>(operator[](r_index)); }

    VectorStatic<T, Row> col(size_t c_index) const {
        VectorStatic<T, Row> result;
        for (size_t i = 0; i < Row; ++i) result[i] = vector_[i][c_index];
        return result;
    }

    VectorStatic<T, Order> diag() const {
        VectorStatic<T, Order> result;
        for (size_t i = 0; i < Order; ++i) result[i] = vector_[i][i];
        return result;
    }

    MatrixStatic<T, Col, Row> transpose() const {
        MatrixStatic<T, Col, Row> result;
        for (size_t i = 0; i < Row; ++i)
            for (size_t j = 0; j < Col; ++j)
                result[j][i] = vector_[i][j];
        return result;
    }

    T det() const {
        static_assert(is_square_(), "Can only find determinant of a square matrix.");
        MatrixStaticLU<T, Order> lu = vt::move(LU());
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
        static_assert(is_square_(), "Can only find trace of a square matrix.");
        T acc = 0;
        for (size_t i = 0; i < Order; ++i) acc += vector_[i][i];
        return acc;
    }

    MatrixStatic inv() const {
        static_assert(is_square_(), "Can only find inverse of a square matrix.");
        MatrixStaticLU<T, Order> lu = vt::move(LU());
        return inv_ut(lu.u()) * inv_lt(lu.l());
    }

    constexpr MatrixStatic inverse() const { return inv(); }

    MatrixStaticLU<T, Order> LU() const {
        static_assert(is_square_(), "Can only find LU decomposition of a square matrix.");
        MatrixStatic<T, Order, Order> lower;
        MatrixStatic<T, Order, Order> upper;
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

    MatrixStatic RRE() const {
        MatrixStatic m(*this);
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
    bool operator==(const MatrixStatic<T, ORow, OCol> &other) const {
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
    bool operator!=(const MatrixStatic<T, ORow, OCol> &other) const { return !operator==(other); }

    template<size_t ORow, size_t OCol>
    bool operator!=(const T (&array)[ORow][OCol]) const { return !operator==(array); }

    template<size_t ORow, size_t OCol>
    bool equals(const MatrixStatic<T, ORow, OCol> &other) const { return operator==(other); }

    template<size_t ORow, size_t OCol>
    bool float_equals(const MatrixStatic<T, ORow, OCol> &other) const {
        if (this == &other) return true;
        if (Row != ORow || Col != OCol) return false;
        for (size_t i = 0; i < Row; ++i) if (!vector_[i].float_equals(other.vector_[i])) return false;
        return true;
    }

    template<size_t ORow, size_t OCol>
    bool equals(const T (&array)[ORow][OCol]) const { return operator==(array); }

    vt::Iterator<VectorStatic<T, Col>> begin() { return vector_.begin(); }

    vt::Iterator<VectorStatic<T, Col>> begin() const { return vector_.begin(); }

    vt::Iterator<VectorStatic<T, Col>> end() { return vector_.end(); }

    vt::Iterator<VectorStatic<T, Col>> end() const { return vector_.end(); }

    constexpr size_t r() const { return Row; }

    constexpr size_t c() const { return Col; }

    constexpr size_t n() const { return Row * Col; }

    constexpr size_t order() const { return Order; }

    constexpr bool is_square() const { return Row == Col; }

private:
    void fix_zero() {
        for (size_t i = 0; i < Row; ++i)
            for (size_t j = 0; j < Col; ++j)
                if (vector_[i][j] == 0.0)
                    vector_[i][j] = 0.0;
    }

    void put_array(size_t index, const VectorStatic<T, Col> &vector) { vector_[index] = vector; }

    void put_array(size_t index, const T (&array)[Col]) { vector_[index] = array; }

    template<typename... Vectors>
    void put_array(size_t index, const VectorStatic<T, Col> &vector, const Vectors &...vectors) {
        vector_[index] = vector;
        if (index < Row - 1) put_array(index + 1, vectors...);
    }

    template<typename... Arrays>
    void put_array(size_t index, const T (&array)[Col], const Arrays (&...arrays)[Col]) {
        vector_[index] = array;
        if (index < Row - 1) put_array(index + 1, arrays...);
    }

    void allocate_zero() { vector_ = vt::move(VectorStatic<VectorStatic<T, Col>, Row>()); }

    void allocate_fill(T fill) { vector_ = vt::move(VectorStatic<VectorStatic<T, Col>, Row>(VectorStatic<T, Col>(fill))); }

    void allocate_from(const MatrixStatic &other) { vector_ = other.vector_; }

    void allocate_from(const T (&array)[Row][Col]) {
        for (size_t i = 0; i < Row; ++i)
            for (size_t j = 0; j < Col; ++j)
                vector_[i][j] = array[i][j];
    }

    void allocate_from(const VectorStatic<T, Col> (&vectors)[Row]) {
        for (size_t i = 0; i < Row; ++i)
            for (size_t j = 0; j < Col; ++j)
                vector_[i][j] = vectors[i][j];
    }

    void steal(MatrixStatic &&other) { vector_ = vt::move(other.vector_); }

    static constexpr bool is_square_() { return Row == Col; }

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
    static MatrixStatic<T, ORow, OCol> &mm_naive(MatrixStatic<T, ORow, OCol> &C,
                                                 const MatrixStatic<T, ORow, X> &A,
                                                 const MatrixStatic<T, X, OCol> &B) {
        for (size_t i = 0; i < ORow; ++i) {
            VectorStatic<T, OCol> &row_C = C[i];
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
    static MatrixStatic<T, ORow, OCol> &mm_strassen(MatrixStatic<T, ORow, OCol> &C,
                                                    const MatrixStatic<T, ORow, X> &A,
                                                    const MatrixStatic<T, X, OCol> &B) {
        // to be implemented later.
        // not suitable for embedded system.
        return mm_naive(C, A, B);
    }

    template<size_t OSize>
    static MatrixStatic<T, OSize, OSize> inv_lt(MatrixStatic<T, OSize, OSize> &L) {
        static_assert(is_square_(), "Can only inverse a square matrix.");
        MatrixStatic<T, OSize, OSize> X;
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
    static MatrixStatic<T, OSize, OSize> inv_ut(MatrixStatic<T, OSize, OSize> &U) {
        static_assert(is_square_(), "Can only inverse a square matrix.");
        MatrixStatic<T, OSize, OSize> tmp = vt::move(U.transpose());
        return inv_lt(tmp).transpose();
    }

public:
    static MatrixStatic zeros() { return MatrixStatic(); }

    static MatrixStatic ones() { return MatrixStatic(1); }

    static MatrixStatic identity() {
        static_assert(is_square_(), "Identity matrix must be a square matrix.");
        return diagonal(1);
    }

    static MatrixStatic id() { return identity(); }

    static MatrixStatic diagonal(T value) {
        MatrixStatic tmp;
        for (size_t i = 0; i < Order; ++i) tmp[i][i] = value;
        return tmp;
    }

    template<size_t ORow, size_t OCol>
    static MatrixStatic quad(const MatrixStatic<T, ORow, OCol> M) { return MatrixStatic(M, M, M, M); }
};

template<typename T, size_t Row, size_t Col>
MatrixStatic<T, Row, Col> operator*(T lhs, const MatrixStatic<T, Row, Col> &rhs) {
    MatrixStatic<T, Row, Col> tmp(rhs);
    tmp.operator*=(lhs);
    return tmp;
}

template<typename T, size_t Row, size_t Col>
T det(const MatrixStatic<T, Row, Col> &A) { return A.det(); }

template<typename T, size_t Row, size_t Col>
T tr(const MatrixStatic<T, Row, Col> &A) { return A.tr(); }

template<typename T, size_t Row, size_t Col>
MatrixStatic<T, Row, Col> inv(const MatrixStatic<T, Row, Col> &A) { return A.inv(); }

template<typename T, size_t Row, size_t Col>
MatrixStatic<T, Row, Col> RRE(const MatrixStatic<T, Row, Col> &A) { return A.RRE(); }

template<typename T, size_t OSize>
class MatrixStaticLU {
private:
    using Matrix_t = MatrixStatic<T, OSize, OSize>;
    using Pair_t = vt::pair<MatrixStatic<T, OSize, OSize>, MatrixStatic<T, OSize, OSize>>;
    MatrixStatic<T, OSize, OSize> l_;
    MatrixStatic<T, OSize, OSize> u_;

public:
    MatrixStaticLU(const MatrixStaticLU &other) : l_(other.l_), u_(other.u_) {}

    MatrixStaticLU(MatrixStaticLU &&other) noexcept: l_(vt::move(other.l_)), u_(vt::move(other.u_)) {}

    explicit MatrixStaticLU(const Pair_t &lu) : l_(lu.first), u_(lu.second) {}

    MatrixStaticLU(const Matrix_t &l, const Matrix_t &u) : l_(l), u_(u) {}

    Matrix_t &l() { return l_; }

    constexpr const Matrix_t &l() const { return l_; }

    Matrix_t &u() { return u_; }

    constexpr const Matrix_t &u() const { return u_; }
};

template<size_t Row, size_t Col = Row>
using numeric_matrix = MatrixStatic<double, Row, Col>;

#endif //VNET_LINALG_STATIC_NUMERIC_MATRIX_H
