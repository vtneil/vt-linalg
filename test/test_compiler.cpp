#include "standard_utility.h"
#pragma GCC optimize (O3)

int x;

void set1() {
    x = 1;
}

void set2() {
    x = 1;
}


int main() {
    vt::if_constexpr<true>::run(set1);
//    vt::if_constexpr<false>::run(set2);

    return x;
}