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
    numeric_vector v({1, 1});

    assert(v.equals({1, 1}));
    assert(v.equals(v));
    assert(numeric_vector().equals({}));
    assert(numeric_vector(5).equals({0, 0, 0, 0, 0}));
    assert(numeric_vector(5, 9).equals({9, 9, 9, 9, 9}));
    assert(numeric_vector(v).equals({1, 1}));
    assert((numeric_vector({1, 2, 3, 4}) + numeric_vector({2, 3, -4, 5})).equals({3, 5, -1, 9}));
    assert((numeric_vector({1, 2, 3, 4}).add({2, 3, -4, 5})).equals({3, 5, -1, 9}));
    assert((numeric_vector({1, 2, 3, 4}) - numeric_vector({2, 3, -4, 5})).equals({-1, -1, 7, -1}));
    assert((numeric_vector({1, 2, 3, 4}).subtract({2, 3, -4, 5})).equals({-1, -1, 7, -1}));
    assert(numeric_vector({1, 2, 3}).dot(numeric_vector({3, 1, -1})) == 2);
    assert(numeric_vector({1, 2, 3}).dot({3, 1, -1}) == 2);

    print(v);

    return 0;
}
