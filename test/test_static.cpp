#include <iostream>
#include <vt_linalg>

template<size_t Size>
using numeric_vector = vt::numeric_vector<Size>;

template<size_t Row, size_t Col = Row>
using numeric_matrix = vt::numeric_matrix<Row, Col>;

#define NEWLINE() std::cout<<'\n'

int main() {
    constexpr size_t N = 128;
    numeric_vector<2> v({1, 1});
    numeric_matrix<3, 4> A({{1, 2, 3, 4},
                            {4, 5, 6, 1},
                            {7, 8, 9, 0}});
    numeric_matrix<4, 2> B({{3, 4},
                            {1, -1},
                            {6, 1},
                            {4, 5}});
    numeric_matrix<3> C({{2, 0, 2},
                         {0, 4, 2},
                         {2, 2, 2}});
    numeric_matrix<N> X;
    numeric_matrix<N> Y;
    numeric_matrix<N> Z1 = vt::move(X * Y);
    numeric_matrix<N> Z2 = vt::move(X.matmul_naive(Y));

    auto v1 = vt::make_numeric_vector({1, 2, 3, 4});
    auto M1 = vt::make_numeric_matrix({{1, 2, 3}, {4, 5, 6}});

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            X[i][j] = static_cast<vt::real_t>(i * N + j) / 762.;
            Y[i][j] = static_cast<vt::real_t>(i * N + j) / 1983.;
        }
    }

    assert(v.equals({1, 1}));
    assert(v.equals(v));
    assert(numeric_vector<5>().equals({0, 0, 0, 0, 0}));
    assert(numeric_vector<5>(9).equals({9, 9, 9, 9, 9}));
    assert(numeric_vector<v.size()>(v).equals({ 1, 1 }));
    assert(numeric_vector<4>({1, 2, 3, 4})[0] == 1);
    assert(numeric_vector<4>({1, 2, 3, 4})[1] == 2);
    assert(numeric_vector<4>({1, 2, 3, 4})[2] == 3);
    assert(numeric_vector<4>({1, 2, 3, 4})[3] == 4);
    assert(numeric_vector<4>({1, 2, 3, 4}).at(3) == 4);
    assert(numeric_vector<4>({1, 2, 3, 4})(3) == 4);
    assert((numeric_vector<4>({1, 2, 3, 4}) + numeric_vector<4>({2, 3, -4, 5})).equals({3, 5, -1, 9}));
    assert((numeric_vector<4>({1, 2, 3, 4}).add({2, 3, -4, 5})).equals({3, 5, -1, 9}));
    assert((numeric_vector<4>({1, 2, 3, 4}) - numeric_vector<4>({2, 3, -4, 5})).equals({-1, -1, 7, -1}));
    assert((numeric_vector<4>({1, 2, 3, 4}).subtract({2, 3, -4, 5})).equals({-1, -1, 7, -1}));
    assert(numeric_vector<3>({1, 2, 3}).inner(numeric_vector<3>({3, 1, -1})) == 2);
    assert(numeric_vector<3>({1, 2, 3}).inner({3, 1, -1}) == 2);
    assert(numeric_vector<3>({1, 2, 3}).outer({3, 2, 1}).equals({{3, 2, 1},
                                                                 {6, 4, 2},
                                                                 {9, 6, 3}}));
    assert(numeric_vector<3>({2, 2, -1}).norm() == 3.);
    assert(numeric_vector<3>({-2, 2, 1}).normalize().equals({-2. / 3., 2. / 3., 1. / 3.}));
    assert(A.equals({{1, 2, 3, 4},
                     {4, 5, 6, 1},
                     {7, 8, 9, 0}}));
    assert((numeric_matrix<A.r(), A.c()>(A) == A));
    assert((numeric_matrix<A.r(), A.c()>(A).equals({{1, 2, 3, 4},
                                                    {4, 5, 6, 1},
                                                    {7, 8, 9, 0}})));
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
    assert((C ^ 9).equals(512. * numeric_matrix<3>({{1160, 2561, 2070},
                                                    {2561, 5791, 4631},
                                                    {2070, 4631, 3721}})));
    assert((C ^ 9).equals(C.matpow(9)));
    assert((C ^ 8).equals(C.matpow_naive(8)));
    assert((C ^ 9).equals(C.matpow_naive(9)));
    assert((C ^ 10).equals(C.matpow_naive(10)));
    assert((C ^ 0).equals(C.matpow_naive(0)));
    assert(A.row(1).equals({4, 5, 6, 1}));
    assert(A.col(2).equals({3, 6, 9}));
    assert(A.diag().equals({1, 5, 9}));
    assert(A.transpose().equals({{1, 4, 7},
                                 {2, 5, 8},
                                 {3, 6, 9},
                                 {4, 1, 0}}));
    assert(det(C) == -8.);
    assert(tr(C) == 8.);
    assert(inv(C) == 0.5 * numeric_matrix<3>({{-1, -1, 2},
                                              {-1, 0,  1},
                                              {2,  1,  -2}}));
    assert(RRE(C) == numeric_matrix<3>::id());
    assert(B * v == numeric_vector<4>({7, 0, 7, 9}));
    assert(Z1.float_equals(Z2));
    assert(v.end() - v.begin() == 2);
    assert((v1.slice<0, 1>().equals({1})));
    assert((v1.slice<0, 4>().equals({1, 2, 3, 4})));
    assert((v1.slice<1, 2>().equals({2})));

    std::cout << "Test Passed!";
    return 0;
}
