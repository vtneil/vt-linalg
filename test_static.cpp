#include <iostream>
#include "vnet_linalg/static_numeric_vector.h"
#include "vnet_linalg/static_numeric_matrix.h"

template<size_t Row, size_t Col = Row>
using num_mat = MatrixStatic<double, Row, Col>;

int main() {
    num_mat<5, 5> A = num_mat<5, 5>::zeros();

    return 0;
}
