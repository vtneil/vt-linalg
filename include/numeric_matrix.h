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
    namespace impl {
        template<typename T, size_t OSize>
        class numeric_matrix_static_lu_t;

        /**
         * Numeric matrix template class where the dimension must be known at compile-time
         * and can't be changed by any ways during runtime to prevent unexpected
         * behavior or errors. This avoid heap usages for safety-critical and
         * performance-critical systems.
         *
         * @tparam T data type
         * @tparam Row row dimension
         * @tparam Col column dimension
         */
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
//        numeric_vector_static_t<T, Col> vector_[Row] = {};

        public:
            /**
             * Default constructor, initializes to zero
             */
            constexpr numeric_matrix_static_t() = default;

            /**
             * Fill constructor, initializes to fill value
             *
             * @param fill Fill value
             */
            constexpr explicit numeric_matrix_static_t(const T &fill)
                    : vector_(
                    numeric_vector_static_t<numeric_vector_static_t<T, Col>, Row>(numeric_vector_static_t<T, Col>(fill))
            ) {}

            /**
             * Copy constructor
             *
             * @param other Other matrix
             */
            constexpr numeric_matrix_static_t(const numeric_matrix_static_t &other) = default;

            /**
             * Move constructor
             *
             * @param other Other vector
             */
            constexpr numeric_matrix_static_t(numeric_matrix_static_t &&other) noexcept = default;

            /**
             * Array constructor, construct from array
             *
             * @param array Array of entries
             */
            constexpr explicit numeric_matrix_static_t(const T (&array)[Row][Col])
                    : vector_(vt::detail::make_nested(array)) {}

            /**
             * Array of vectors as rows constructor, construct from array
             *
             * @param vectors Array of rows
             */
            explicit numeric_matrix_static_t(const numeric_vector_static_t<T, Col> (&vectors)[Row]) {
                allocate_from(vectors);
            }

            /**
             * Nested vector constructor
             *
             * @param nested Nested vector
             */
            constexpr explicit numeric_matrix_static_t(
                    const numeric_vector_static_t<numeric_vector_static_t<T, Col>, Row> &nested)
                    : vector_(nested) {}

            /**
             * Block matrix constructor
             *
             * @tparam ORow
             * @tparam OCol
             * @tparam M
             * @tparam N
             * @param blocks Array of blocks
             */
            template<size_t ORow, size_t OCol, size_t M, size_t N>
            explicit numeric_matrix_static_t(const numeric_matrix_static_t<T, ORow, OCol> (&blocks)[M][N]) {
                helper_insert_major(0, blocks);
            }

            /**
             * Quad-matrix constructor, construct M from M11, M12, M21, M22
             *
             * @tparam R1
             * @tparam C1
             * @param M11 Upper-left matrix
             * @param M12 Upper-right matrix
             * @param M21 Lower-left matrix
             * @param M22 Lower-right matrix
             */
            template<size_t R1, size_t C1>
            numeric_matrix_static_t(const numeric_matrix_static_t<T, R1, C1> &M11,
                                    const numeric_matrix_static_t<T, R1, Col - C1> &M12,
                                    const numeric_matrix_static_t<T, Row - R1, C1> &M21,
                                    const numeric_matrix_static_t<T, Row - R1, Col - C1> &M22) {
                insert<0, 0>(M11);
                insert<0, C1>(M12);
                insert<R1, 0>(M21);
                insert<R1, C1>(M22);
            }

            /**
             * B matrix constructor, construct M from [A | B]
             *
             * @tparam C1 Main column
             * @param A Main matrix
             * @param B Augment part
             */
            template<size_t C1>
            numeric_matrix_static_t(const numeric_matrix_static_t<T, Row, C1> &A,
                                    const numeric_matrix_static_t<T, Row, Col - C1> &B) {
                insert<0, 0>(A);
                insert<0, C1>(B);
            }

            FORCE_INLINE numeric_vector_static_t<T, Col> &operator[](size_t index) { return vector_[index]; }

            FORCE_INLINE constexpr const numeric_vector_static_t<T, Col> &
            operator[](size_t index) const { return vector_[index]; }

            FORCE_INLINE T &at(size_t r_index, size_t c_index) { return vector_[r_index][c_index]; };

            FORCE_INLINE constexpr const T &
            at(size_t r_index, size_t c_index) const { return vector_[r_index][c_index]; };

            FORCE_INLINE T &operator()(size_t r_index, size_t c_index) { return at(r_index, c_index); }

            FORCE_INLINE constexpr const T &operator()(size_t r_index, size_t c_index) const {
                return at(r_index, c_index);
            }

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

            numeric_matrix_static_t &operator+=(const numeric_matrix_static_t &other) {
                return iadd(*this, *this, other);
            }

            numeric_matrix_static_t operator+(const numeric_matrix_static_t &other) const {
                numeric_matrix_static_t tmp(*this);
                tmp.operator+=(other);
                return tmp;
            }

            constexpr numeric_matrix_static_t add(const numeric_matrix_static_t &other) const {
                return operator+(other);
            }

            /**
             * In-place C = A + B_ tools
             *
             * @param C
             * @param A
             * @param B
             * @return C
             */
            static numeric_matrix_static_t &iadd(numeric_matrix_static_t &C,
                                                 const numeric_matrix_static_t &A,
                                                 const numeric_matrix_static_t &B) {
                for (size_t i = 0; i < Row; ++i)
                    for (size_t j = 0; j < Col; ++j)
                        C[i][j] = A[i][j] + B[i][j];
                return C;
            }

            numeric_matrix_static_t &operator-=(const numeric_matrix_static_t &other) {
                return isub(*this, *this, other);
            }

            numeric_matrix_static_t operator-(const numeric_matrix_static_t &other) const {
                numeric_matrix_static_t tmp(*this);
                tmp.operator-=(other);
                return tmp;
            }

            constexpr numeric_matrix_static_t sub(const numeric_matrix_static_t &other) const {
                return operator-(other);
            }

            /**
             * In-place C = A - B_ tools
             *
             * @param C
             * @param A
             * @param B
             * @return C
             */
            static numeric_matrix_static_t &isub(numeric_matrix_static_t &C,
                                                 const numeric_matrix_static_t &A,
                                                 const numeric_matrix_static_t &B) {
                for (size_t i = 0; i < Row; ++i)
                    for (size_t j = 0; j < Col; ++j)
                        C[i][j] = A[i][j] - B[i][j];
                return C;
            }

            numeric_matrix_static_t<T, Order, Order> &
            operator*=(const numeric_matrix_static_t<T, Order, Order> &other) {
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

            /**
             * In-place C = AB tools
             *
             * @tparam ORow
             * @tparam X
             * @tparam OCol
             * @param C
             * @param A
             * @param B
             * @return C
             */
            template<size_t ORow, size_t X, size_t OCol>
            static constexpr numeric_matrix_static_t<T, ORow, OCol> imatmul(numeric_matrix_static_t<T, ORow, OCol> &C,
                                                                            const numeric_matrix_static_t<T, ORow, X> &A,
                                                                            const numeric_matrix_static_t<T, X, OCol> &B) {
                return mm_naive(C, A, B);
            }

            /**
             * Multiply this matrix with the other matrix in naive O(n^3) way (default multiplication method).
             *
             * @tparam ORow
             * @tparam OCol
             * @param other
             * @return Multiplied matrix
             */
            template<size_t ORow, size_t OCol>
            numeric_matrix_static_t<T, Row, OCol>
            matmul_naive(const numeric_matrix_static_t<T, ORow, OCol> &other) const {
                numeric_matrix_static_t<T, Row, OCol> C;
                return mm_naive(C, *this, other);
            }

            /**
             * Transform input vector to output vector by left-multiplication of this matrix B_ = Ax.\n
             * Alias for writing A * x_ (which is more preferred than transform function.)
             *
             * @param other
             * @return Transformed vector
             */
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
                if (n == 0) return identity();
                if (n == 1) return numeric_matrix_static_t(*this);
                numeric_matrix_static_t<T, Order, Order> base(*this);
                numeric_matrix_static_t<T, Order, Order> product(vt::move(identity()));
                while (n > 0) {
                    if (n % 2 == 1) product.operator*=(base);
                    if (n > 1) base *= base;
                    n >>= 1;
                }
                return product;
            }

            /**
             * Calculate a nth power of this matrix in O(log n * mulO). Must be a square matrix.\n
             * If this matrix is not square, the compile-time error is thrown.
             *
             * @param n
             * @return A^n
             */
            numeric_matrix_static_t<T, Order, Order> matpow(size_t n) {
                static_assert(static_is_a_square_matrix(), "Non-square matrix can\'t use power operator.");
                return operator^(n);
            }

            /**
             * Calculate a nth power of this matrix in O(n * mulO) naively. Must be a square matrix.\n
             * If this matrix is not square, the compile-time error is thrown.
             *
             * @param n
             * @return A^n
             */
            numeric_matrix_static_t<T, Order, Order> matpow_naive(size_t n) {
                static_assert(static_is_a_square_matrix(), "Non-square matrix can\'t use power operator.");
                numeric_matrix_static_t<T, Order, Order> product = vt::move(id());
                for (size_t i = 0; i < n; ++i) product *= *this;
                return product;
            }

            /**
             * Insert the other matrix into current matrix, returning reference to the current matrix.\n
             * If insertion is not valid, the compile-time error is thrown.
             *
             * @tparam pos_row Row position to insert, default to 0
             * @tparam pos_col Column position to insert, default to 0
             * @tparam ORow Other matrix's row size
             * @tparam OCol Other matrix's column size
             * @param M Other matrix
             * @return Reference to this matrix
             */
            template<size_t pos_row = 0, size_t pos_col = 0, size_t ORow, size_t OCol>
            numeric_matrix_static_t &insert(const numeric_matrix_static_t<T, ORow, OCol> &M) {
                static_assert(pos_row < Row, "Insertion failed! Start row is out of range.");
                static_assert(pos_col < Col, "Insertion failed! Start column is out of range.");
                static_assert(pos_row + ORow <= Row, "Insertion failed! Stop row is out of range.");
                static_assert(pos_col + OCol <= Col, "Insertion failed! Stop column is out of range.");
                for (size_t i = 0; i < ORow; ++i)
                    for (size_t j = 0; j < OCol; ++j)
                        vector_[pos_row + i][pos_col + j] = M[i][j];
                return *this;
            }

            /**
             * Insert the other matrix into current matrix, returning reference to the current matrix.\n
             * If insertion is not valid, the compile-time error is thrown.
             *
             * @tparam pos_row Row position to insert, default to 0
             * @tparam pos_col Column position to insert, default to 0
             * @tparam ORow Other matrix's row size
             * @tparam OCol Other matrix's column size
             * @param array Other matrix as 2D array
             * @return Reference to this matrix
             */
            template<size_t pos_row = 0, size_t pos_col = 0, size_t ORow, size_t OCol>
            numeric_matrix_static_t &insert(const T (&array)[ORow][OCol]) {
                static_assert(pos_row < Row, "Insertion failed! Start row is out of range.");
                static_assert(pos_col < Col, "Insertion failed! Start column is out of range.");
                static_assert(pos_row + ORow <= Row, "Insertion failed! Stop row is out of range.");
                static_assert(pos_col + OCol <= Col, "Insertion failed! Stop column is out of range.");
                for (size_t i = 0; i < ORow; ++i)
                    for (size_t j = 0; j < OCol; ++j)
                        vector_[pos_row + i][pos_col + j] = array[i][j];
                return *this;
            }

            /**
             * Slice current matrix.
             *
             * @tparam r1 Row position from
             * @tparam c1 Column position from
             * @tparam r2 Row position to
             * @tparam c2 Colum position to
             * @return A slice of current matrix
             */
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

            /**
             * Returns row at index as vector.\n
             * WARNING: Doesn't have out of range check!
             *
             * @param r_index Row index
             * @return Row at index as vector
             */
            constexpr numeric_vector_static_t<T, Col> row(size_t r_index) const {
                return numeric_vector_static_t<T, Col>(operator[](r_index));
            }

            /**
             * Returns column at index as vector.\n
             * WARNING: Doesn't have out of range check!
             *
             * @param c_index Column index
             * @return Column at index as vector
             */
            numeric_vector_static_t<T, Row> col(size_t c_index) const {
                numeric_vector_static_t<T, Row> result;
                for (size_t i = 0; i < Row; ++i) result[i] = vector_[i][c_index];
                return result;
            }

            /**
             * Returns main diagonal entries as vector.
             *
             * @return Main diagonal entries as vector
             */
            numeric_vector_static_t<T, Order> diag() const {
                numeric_vector_static_t<T, Order> result;
                for (size_t i = 0; i < Order; ++i) result[i] = vector_[i][i];
                return result;
            }

            /**
             * Returns A^T (transposed matrix of this matrix).
             *
             * @return A^T
             */
            numeric_matrix_static_t<T, Col, Row> transpose() const {
                numeric_matrix_static_t<T, Col, Row> result;
                for (size_t i = 0; i < Row; ++i)
                    for (size_t j = 0; j < Col; ++j)
                        result[j][i] = vector_[i][j];
                return result;
            }

            /**
             * Finds determinant of this matrix.\n
             * If this matrix is not square, the compile-time error is thrown.
             *
             * @return Determinant of this matrix
             */
            T det() const {
                static_assert(static_is_a_square_matrix(), "Can only find determinant of a square matrix.");
                return det_from_lu(LU());
            }

        private:
            T det_from_lu(const numeric_matrix_static_lu_t<T, Order> &lu) const {
                T det = 1;
                size_t r_swaps = 0;
                for (size_t i = 0; i < Order; ++i) {
                    det *= lu.u()[i][i];
                    if (lu.u()[i][i] == 0) return 0;
                    if (lu.l()[i][i] != 1) r_swaps++;
                }
                return (r_swaps % 2 == 0) ? det : -det;
            }

        public:

            /**
             * Finds trace of this matrix.\n
             * If this matrix is not square, the compile-time error is thrown.
             *
             * @return Trace of this matrix
             */
            T tr() const {
                static_assert(static_is_a_square_matrix(), "Can only find trace of a square matrix.");
                T acc = 0;
                for (size_t i = 0; i < Order; ++i) acc += vector_[i][i];
                return acc;
            }

            /**
             * Finds inverse of this matrix.\n
             * If this matrix is not square, the compile-time error is thrown.
             *
             * If inverse doesn't exist (det = 0), the zero matrix is returned.
             *
             * @return Inverse of this matrix
             */
            numeric_matrix_static_t inv() const {
                static_assert(static_is_a_square_matrix(), "Can only find inverse of a square matrix.");
                numeric_matrix_static_lu_t<T, Order> lu = vt::move(LU());
                if (abs(det_from_lu(lu)) > 1e-10) return (inv_ut(lu.u()) * inv_lt(lu.l()));
                return numeric_matrix_static_t();
            }

            /**
             * Finds inverse of this matrix.\n
             * If this matrix is not square, the compile-time error is thrown.
             *
             * @return Inverse of this matrix
             */
            constexpr numeric_matrix_static_t inverse() const { return inv(); }

            /**
             * Finds LU-decomposition of this matrix.\n
             * If this matrix is not square, the compile-time error is thrown.
             *
             * @return LU-decomposition of this matrix
             */
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

            /**
             * Finds Row-Reduced Echlon (RRE) form of this matrix.\n
             * If this matrix is not square, the compile-time error is thrown.
             *
             * @return Row-Reduced Echlon form of this matrix
             */
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

            /**
             * Checks equality of this matrix and the other matrix.
             *
             * @tparam ORow
             * @tparam OCol
             * @param other Other matrix
             * @return
             */
            template<size_t ORow, size_t OCol>
            bool operator==(const numeric_matrix_static_t<T, ORow, OCol> &other) const {
                if (this == &other) return true;
                if (Row != ORow || Col != OCol) return false;
                for (size_t i = 0; i < Row; ++i) if (vector_[i] != other.vector_[i]) return false;
                return true;
            }

            /**
             * Checks equality of this matrix and the other matrix.
             *
             * @tparam ORow
             * @tparam OCol
             * @param array Other matrix as array
             * @return
             */
            template<size_t ORow, size_t OCol>
            bool operator==(const T (&array)[ORow][OCol]) const {
                if (Row != ORow || Col != OCol) return false;
                for (size_t i = 0; i < Row; ++i) if (vector_[i] != array[i]) return false;
                return true;
            }

            /**
             * Checks inequality of this matrix and the other matrix.
             *
             * @tparam ORow
             * @tparam OCol
             * @param other Other matrix
             * @return
             */
            template<size_t ORow, size_t OCol>
            constexpr bool operator!=(const numeric_matrix_static_t<T, ORow, OCol> &other) const {
                return !operator==(other);
            }

            /**
             * Checks inequality of this matrix and the other matrix.
             *
             * @tparam ORow
             * @tparam OCol
             * @param array Other matrix as array
             * @return
             */
            template<size_t ORow, size_t OCol>
            constexpr bool operator!=(const T (&array)[ORow][OCol]) const { return !operator==(array); }

            /**
             * Checks equality of this matrix and the other matrix.
             *
             * @tparam ORow
             * @tparam OCol
             * @param other Other matrix
             * @return
             */
            template<size_t ORow, size_t OCol>
            constexpr bool equals(const numeric_matrix_static_t<T, ORow, OCol> &other) const {
                return operator==(other);
            }

            /**
             * Checks equality of this matrix and the other matrix with float/double threshold using
             * equation abs(x_ - y) < threshold for equality.
             * @tparam ORow
             * @tparam OCol
             * @param other Other matrix
             * @return
             */
            template<size_t ORow, size_t OCol>
            bool float_equals(const numeric_matrix_static_t<T, ORow, OCol> &other) const {
                if (this == &other) return true;
                if (Row != ORow || Col != OCol) return false;
                for (size_t i = 0; i < Row; ++i) if (!vector_[i].float_equals(other.vector_[i])) return false;
                return true;
            }

            /**
             * Checks equality of this matrix and the other matrix.
             *
             * @tparam ORow
             * @tparam OCol
             * @param array Other matrix as array
             * @return
             */
            template<size_t ORow, size_t OCol>
            constexpr bool equals(const T (&array)[ORow][OCol]) const { return operator==(array); }

            /**
             * Returns an iterator to matrix's first row.
             *
             * @return An iterator to matrix's first row
             */
            iterator<numeric_vector_static_t<T, Col>> begin() { return vector_.begin(); }

            /**
             * Returns an iterator to matrix's first row.
             *
             * @return An iterator to matrix's first row
             */
            constexpr iterator<numeric_vector_static_t<T, Col>> begin() const { return vector_.begin(); }

            /**
             * Returns an iterator to matrix's last row.
             *
             * @return An iterator to matrix's last row
             */
            iterator<numeric_vector_static_t<T, Col>> end() { return vector_.end(); }

            /**
             * Returns an iterator to matrix's last row.
             *
             * @return An iterator to matrix's last row
             */
            constexpr iterator<numeric_vector_static_t<T, Col>> end() const { return vector_.end(); }

            /**
             * Returns number of rows.
             *
             * @return Number of rows
             */
            constexpr size_t r() const { return Row; }

            /**
             * Returns number of columns.
             *
             * @return Number of columns
             */
            constexpr size_t c() const { return Col; }

            /**
             * Returns number of entries.
             *
             * @return Number of entries
             */
            constexpr size_t n() const { return Row * Col; }

            /**
             * Returns order of this matrix (min(Row, Col)).
             *
             * @return Order of this matrix
             */
            constexpr size_t order() const { return Order; }

            /**
             * Returns whether the matrix is square.
             *
             * @return
             */
            constexpr bool is_square() const { return Row == Col; }

            /**
             * Swaps entries with the other matrix.
             *
             * @param other
             */
            void swap(numeric_matrix_static_t &other) {
                for (size_t i = 0; i < Row; ++i) vector_[i].swap(other.vector_[i]);
            }

            /**
             * Creates a copy of this matrix.
             *
             * @return A copy of this matrix
             */
            constexpr numeric_matrix_static_t copy() const { return numeric_matrix_static_t(*this); }

        private:
            void fix_zero() {
                for (size_t i = 0; i < Row; ++i)
                    for (size_t j = 0; j < Col; ++j)
                        if (vector_[i][j] == 0.0)
                            vector_[i][j] = 0.0;
            }

            void allocate_zero() {
                vector_ = vt::move(numeric_vector_static_t<numeric_vector_static_t<T, Col>, Row>());
            }

            void allocate_fill(T fill) {
                vector_ = vt::move(
                        numeric_vector_static_t<numeric_vector_static_t<T, Col>, Row>(
                                numeric_vector_static_t<T, Col>(fill))
                );
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

            template<size_t ORow, size_t OCol, size_t M, size_t N>
            void helper_insert_major(size_t pos_row, const numeric_matrix_static_t<T, ORow, OCol> (&blocks)[M][N]) {
                if (pos_row < M) {
                    helper_insert_minor(pos_row, 0, blocks);
                    helper_insert_major(pos_row + 1, blocks);
                }
            }

            template<size_t ORow, size_t OCol, size_t M, size_t N>
            void helper_insert_minor(size_t pos_row, size_t pos_col,
                                     const numeric_matrix_static_t<T, ORow, OCol> (&blocks)[M][N]) {
                helper_insert_unsafe(ORow * pos_row, OCol * pos_col, blocks[pos_row][pos_col]);
                if (pos_col < N - 1) helper_insert_minor(pos_row, pos_col + 1, blocks);
            }

            template<size_t ORow, size_t OCol>
            void helper_insert_unsafe(size_t pos_row, size_t pos_col, const numeric_matrix_static_t<T, ORow, OCol> &M) {
                for (size_t i = 0; i < ORow; ++i)
                    for (size_t j = 0; j < OCol; ++j)
                        vector_[pos_row + i][pos_col + j] = M[i][j];
            }

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
            static constexpr numeric_matrix_static_t<T, ORow, OCol> &
            mm_strassen(numeric_matrix_static_t<T, ORow, OCol> &C,
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
            /**
             * Returns zero matrix.
             *
             * @return Zero matrix
             */
            static constexpr numeric_matrix_static_t zeros() { return numeric_matrix_static_t(); }

            /**
             * Returns 1-filled matrix.
             *
             * @return 1-filled matrix
             */
            static constexpr numeric_matrix_static_t ones() { return numeric_matrix_static_t(1); }

            /**
             * Returns an identity matrix of current dimension.\n
             * If this matrix is not square, the compile-time error is thrown.
             *
             * @return Identity matrix
             */
            static constexpr numeric_matrix_static_t identity() {
                static_assert(static_is_a_square_matrix(), "Identity matrix must be a square matrix.");
                return diagonals(1);
            }

            /**
             * Returns an identity matrix of current dimension.\n
             * If this matrix is not square, the compile-time error is thrown.
             *
             * @return Identity matrix
             */
            static constexpr numeric_matrix_static_t id() { return identity(); }

            /**
             * Returns a matrix where the main diagonal is filled with fill value.
             *
             * @param value Fill value
             * @return Filled diagonal matrix
             */
            static numeric_matrix_static_t diagonals(T value) {
                numeric_matrix_static_t tmp;
                for (size_t i = 0; i < Order; ++i) tmp[i][i] = value;
                return tmp;
            }

            /**
             * Returns a matrix where the main diagonal is filled with fill values from array.
             *
             * @param array Fill array
             * @return Filled diagonal matrix
             */
            static numeric_matrix_static_t diagonals(const T (&array)[Order]) {
                numeric_matrix_static_t tmp;
                for (size_t i = 0; i < Order; ++i) tmp[i][i] = array[i];
                return tmp;
            }

            /**
             * Returns a quad-filled matrix where M11, M12, M21, and M22 are the same.
             *
             * @tparam ORow
             * @tparam OCol
             * @param M Fill matrix
             * @return Quad-filled matrix
             */
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
    }

    /**
     * Finds determinant of this matrix.\n
     * If this matrix is not square, the compile-time error is thrown.
     *
     * @tparam T
     * @tparam Row
     * @tparam Col
     * @param A
     * @return Determinant
     */
    template<typename T, size_t Row, size_t Col>
    T det(const impl::numeric_matrix_static_t <T, Row, Col> &A) { return A.det(); }

    /**
     * Finds trace of this matrix.\n
     * If this matrix is not square, the compile-time error is thrown.
     *
     * @tparam T
     * @tparam Row
     * @tparam Col
     * @param A
     * @return Trace
     */
    template<typename T, size_t Row, size_t Col>
    T tr(const impl::numeric_matrix_static_t <T, Row, Col> &A) { return A.tr(); }

    /**
     * Finds inverse of this matrix.\n
     * If this matrix is not square, the compile-time error is thrown.
     *
     * @tparam T
     * @tparam Row
     * @tparam Col
     * @param A
     * @return Inverse
     */
    template<typename T, size_t Row, size_t Col>
    impl::numeric_matrix_static_t <T, Row, Col> inv(const impl::numeric_matrix_static_t <T, Row, Col> &A) {
        return A.inv();
    }

    /**
     * Finds RRE form of this matrix.
     *
     * @tparam T
     * @tparam Row
     * @tparam Col
     * @param A
     * @return RRE form of this matrix
     */
    template<typename T, size_t Row, size_t Col>
    impl::numeric_matrix_static_t <T, Row, Col> RRE(const impl::numeric_matrix_static_t <T, Row, Col> &A) {
        return A.RRE();
    }

    namespace impl {
        /**
         * Wrapper class for LU-decomposed matrix comprised of L and U square matrices.
         *
         * @tparam T
         * @tparam OSize
         */
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

            /**
             * L Matrix
             *
             * @return L Matrix
             */
            Matrix_t &l() { return l_; }

            /**
             * L Matrix
             *
             * @return L Matrix
             */
            constexpr const Matrix_t &l() const { return l_; }

            /**
             * U Matrix
             *
             * @return u Matrix
             */
            Matrix_t &u() { return u_; }

            /**
             * U Matrix
             *
             * @return u Matrix
             */
            constexpr const Matrix_t &u() const { return u_; }
        };
    }

    /**
     * Numeric matrix where the dimension must be known at compile-time
     * and can't be changed by any ways during runtime to prevent unexpected
     * behavior or errors. This avoid heap usages for safety-critical and
     * performance-critical systems.
     * \n
     * This numeric matrix uses real type (real_t), defined as double).
     *
     * @tparam T Data type
     * @tparam Row Row dimension
     * @tparam Col Column dimension
     */
    template<size_t Row, size_t Col = Row>
    using numeric_matrix = impl::numeric_matrix_static_t<real_t, Row, Col>;

    template<size_t OSize>
    using numeric_matrix_lu = impl::numeric_matrix_static_lu_t<real_t, OSize>;

    /**
     * Creates a numeric matrix.
     *
     * @tparam Row Row dimension
     * @tparam Col Column dimension
     * @param array Array of data
     * @return Numeric matrix
     */
    template<size_t Row, size_t Col = Row>
    constexpr numeric_matrix<Row, Col> make_numeric_matrix(const real_t (&array)[Row][Col]) {
        return numeric_matrix<Row, Col>(array);
    }

    /**
     * Creates a numeric matrix (copy).
     * @tparam Row Row dimension
     * @tparam Col Column dimension
     * @param M Input numeric matrix
     * @return Numeric matrix
     */
    template<size_t Row, size_t Col = Row>
    constexpr numeric_matrix<Row, Col> make_numeric_matrix(const numeric_matrix<Row, Col> &M) {
        return numeric_matrix<Row, Col>(M);
    }

    /**
     * Creates a quad matrix.
     *
     * @tparam Row
     * @tparam Col
     * @tparam R1
     * @tparam R2
     * @tparam C1
     * @tparam C2
     * @param M11
     * @param M12
     * @param M21
     * @param M22
     * @return Numeric matrix
     */
    template<size_t R1, size_t R2, size_t C1, size_t C2>
    constexpr numeric_matrix<R1 + R2, C1 + C2> make_quad_matrix(const numeric_matrix<R1, C1> &M11,
                                                                const numeric_matrix<R1, C2> &M12,
                                                                const numeric_matrix<R2, C1> &M21,
                                                                const numeric_matrix<R2, C2> &M22) {
        return numeric_matrix<R1 + R2, C1 + C2>(M11, M12, M21, M22);
    }

    /**
     * Creates a block matrix.
     *
     * @tparam ORow
     * @tparam OCol
     * @tparam M
     * @tparam N
     * @param blocks Array of matrices
     * @return Numeric matrix
     */
    template<size_t ORow, size_t OCol, size_t M, size_t N>
    constexpr numeric_matrix<(ORow * M), (OCol * N)>
    make_block_matrix(const numeric_matrix<ORow, OCol> (&blocks)[M][N]) {
        return numeric_matrix<(ORow * M), (OCol * N)>(blocks);
    }
}

#endif //VNET_LINALG_NUMERIC_MATRIX_H
