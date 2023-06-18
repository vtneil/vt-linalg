#include <iostream>
#include "../src/complex_number.h"

using namespace vt;

decltype(std::cout) &print_complex(complex_t number) {
    return std::cout << number.real() << " + " << number.imag() << "i";
}

int main() {
    complex_t a = 4;
    complex_t b(1, 3);

    a *= b;
    print_complex(a) << '\n';

    a = complex_t(3, 4) / complex_t(8, -2);
    a++;

    print_complex(a) << '\n';

    return 0;
}