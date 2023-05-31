#include <iostream>
#include "vnet_linalg/numeric_vector.h"
#include "vnet_linalg/numeric_matrix.h"

#define NEWLINE() std::cout<<'\n'

template<typename T>
void print(const T &x) {
    std::cout << x << '\n';
}

template<typename T>
void print(const Matrix<T> &A) {
    if (A.r() > 0) {
        for (auto &x: A) {
            if (x.size() > 0) {
                for (auto &y: x) {
                    std::cout << y << ' ';
                }
                NEWLINE();
            }
        }
    }
}

template<typename T>
void printb(const Matrix<T> &A) {
    std::cout << '{';
    for (auto &x: A) {
        std::cout << '{';
        for (auto &y: x) {
            std::cout << y;
            if (&y != (x->end() - 1)) std::cout << ',';
        }
        std::cout << '}';
        if (&x != (A.end() - 1)) std::cout << ',';
    }
    std::cout << '}';
    NEWLINE();
}

template<typename T>
void print(const Vector<T> &v) {
    for (auto &y: v) {
        std::cout << y << ' ';
    }
    if (v.size()) NEWLINE();
}

template<typename T>
void print(const Vector<Vector<T>> &v) {
    for (auto &x: v) {
        for (auto &y: x) {
            std::cout << y << ' ';
        }
        if (x.size()) NEWLINE();
    }
    if (v.size()) NEWLINE();
}

int main() {
    constexpr size_t N = 1024;
    numeric_vector v({1, 1});
    numeric_matrix A({{1, 2, 3, 4},
                      {4, 5, 6, 1},
                      {7, 8, 9, 0}});
    numeric_matrix B({{3, 4},
                      {1, -1},
                      {6, 1},
                      {4, 5}});
    numeric_matrix C({{2, 0, 2},
                      {0, 4, 2},
                      {2, 2, 2}});
    numeric_matrix X(N, N);
    numeric_matrix Y(N, N);
    numeric_matrix Z1 = move(X * Y);
    numeric_matrix Z2 = move(X.matmul_naive(Y));

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            X[i][j] = (double) (i * N + j) / 762;
            Y[i][j] = (double) (i * N + j) / 1983;
        }
    }

    assert(v.equals({1, 1}));
    assert(v.equals(v));
    assert(numeric_vector().equals({}));
    assert(numeric_vector(5).equals({0, 0, 0, 0, 0}));
    assert(numeric_vector(5, 9).equals({9, 9, 9, 9, 9}));
    assert(numeric_vector(v).equals({ 1, 1 }));
    assert(numeric_vector({1, 2, 3, 4})[0] == 1);
    assert(numeric_vector({1, 2, 3, 4})[1] == 2);
    assert(numeric_vector({1, 2, 3, 4})[2] == 3);
    assert(numeric_vector({1, 2, 3, 4})[3] == 4);
    assert(numeric_vector({1, 2, 3, 4}).at(3) == 4);
    assert(numeric_vector({1, 2, 3, 4})(3) == 4);
    assert((numeric_vector({1, 2, 3, 4}) + numeric_vector({2, 3, -4, 5})).equals({3, 5, -1, 9}));
    assert((numeric_vector({1, 2, 3, 4}).add({2, 3, -4, 5})).equals({3, 5, -1, 9}));
    assert((numeric_vector({1, 2, 3, 4}) - numeric_vector({2, 3, -4, 5})).equals({-1, -1, 7, -1}));
    assert((numeric_vector({1, 2, 3, 4}).subtract({2, 3, -4, 5})).equals({-1, -1, 7, -1}));
    assert(numeric_vector({1, 2, 3}).inner(numeric_vector({3, 1, -1})) == 2);
    assert(numeric_vector({1, 2, 3}).inner({3, 1, -1}) == 2);
    assert(numeric_vector({1, 2, 3}).outer({3, 2, 1}).equals({{3, 2, 1},
                                                              {6, 4, 2},
                                                              {9, 6, 3}}));
    assert(numeric_vector({2, 2, -1}).norm() == 3.);
    assert(numeric_vector({-2, 2, 1}).normalize().equals({-2. / 3., 2. / 3., 1. / 3.}));
    assert(A.equals({{1, 2, 3, 4},
                     {4, 5, 6, 1},
                     {7, 8, 9, 0}}));
    assert(numeric_matrix(A) == A);
    assert(numeric_matrix(A).equals({{ 1, 2, 3, 4 },
           { 4, 5, 6, 1 },
           { 7, 8, 9, 0 }}));
    assert(A[0][0] == 1);
    assert(A[0][1] == 2);
    assert(A[0][3] == 4);
    assert(A[2][3] == 0);
    assert(A(0, 0) == 1);
    assert(A(0, 1) == 2);
    assert(A(0, 3) == 4);
    assert(A(2, 3) == 0);
    assert(A.at(2, 3) == 0);
    assert(B + B == 2. * B);
    assert((B + B).equals({{6,  8},
                           {2,  -2},
                           {12, 2},
                           {8,  10}}));
    assert(!A.is_square());
    assert(C.is_square());
    assert(A.can_multiply_with(B));
    assert(!A.can_multiply_with(A));
    assert(B.can_multiply_with(v));
    assert(!B.can_multiply_with(numeric_vector({1, 2, 3, 4})));
    assert((C ^ 9).equals(512. * numeric_matrix({{1160, 2561, 2070},
                                                 {2561, 5791, 4631},
                                                 {2070, 4631, 3721}})));
    assert((C ^ 9).equals(C.matpow(9)));
    assert((C ^ 9).equals(C.matpow_naive(9)));
    assert(A.row(1).equals({4, 5, 6, 1}));
    assert(A.col(2).equals({3, 6, 9}));
    assert(A.diag().equals({1, 5, 9}));
    assert(A.transpose().equals({{1, 4, 7},
                                 {2, 5, 8},
                                 {3, 6, 9},
                                 {4, 1, 0}}));
    assert(det(C) == -8.);
    assert(tr(C) == 8.);
    assert(inv(C) == 0.5 * numeric_matrix({{-1, -1, 2},
                                           {-1, 0,  1},
                                           {2,  1,  -2}}));
    assert(RRE(C) == numeric_matrix::id(3));
    assert(B * v == numeric_vector({7, 0, 7, 9}));
    assert(Z1.float_equals(Z2));
    return 0;
}
