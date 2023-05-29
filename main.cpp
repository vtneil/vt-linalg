#include <iostream>
//#include <random>
#include "vnet_linalg/vector.h"
#include "vnet_linalg/matrix.h"

//#define TEST
#define NEWLINE() std::cout<<'\n'

using matrix = Matrix<double>;
using vector = Vector<double>;

template<typename T>
void print(const T &x) {
    std::cout << x << '\n';
}

template<typename T>
void print(const Matrix<T> &A) {
    for (auto &x: A) {
        for (auto &y: *x) {
            std::cout << y << ' ';
        }
        NEWLINE();
    }
}

template<typename T>
void printb(const Matrix<T> &A) {
    std::cout << '{';
    for (auto &x: A) {
        std::cout << '{';
        for (auto &y: *x) {
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

void test(size_t n) {
    size_t N = n;

    matrix A(N, N);
    matrix B(N, N);
    matrix C(N, N);

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            A[i][j] = (double) (i * N + j) / 762;
            B[i][j] = (double) (i * N + j) / 1983;
        }
    }

    for (size_t i = 0; i < 1; ++i) {
        static_cast<void>(matrix::imatmul(C, A, B));
    }
}

int main(int argc, char **argv) {
    matrix A({{1, 2, 3, 4},
              {4, 5, 6, 1},
              {7, 8, 9, 0}});
    matrix B({{3, 4},
              {1, -1},
              {6, 1},
              {4, 5}});

//    A = {{1, 2, 3, 4}, {4, 5, 6, 1}, {7, 8, 9, 0}, {1, 4, -1, 9}};
//    B = {{1, 2, 3, 4}, {4, 5, 6, 1}, {7, 8, 9, 0}, {1, 4, -1, 9}};

    matrix C(A * B);

    print(C);

#ifdef TEST
    if (argc < 2) return -1;
    size_t n = strtoul(argv[1], nullptr, 0);
    test(n);
#endif
    return 0;
}
