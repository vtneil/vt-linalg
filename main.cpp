#include <iostream>
//#include <random>
#include "vnet_linalg/numeric_vector.h"
#include "vnet_linalg/numeric_matrix.h"

//#define TEST
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

void test(size_t n) {
    size_t N = n;

    numeric_matrix A(N, N);
    numeric_matrix B(N, N);
    numeric_matrix C(N, N);

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            A[i][j] = (double) (i * N + j) / 762;
            B[i][j] = (double) (i * N + j) / 1983;
        }
    }

    for (size_t i = 0; i < 1; ++i) {
        static_cast<void>(numeric_matrix::imatmul(C, A, B));
    }

    static_cast<void>(A.r());
}

int main(int argc, char **argv) {
    std::cout << numeric_matrix(5, 8).r();


#ifdef TEST
    if (argc < 2) return -1;
    size_t n = strtoul(argv[1], nullptr, 0);
    test(n);
#endif
    return 0;
}
