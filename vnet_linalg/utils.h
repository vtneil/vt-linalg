#ifndef VNET_LINALG_UTILS_H
#define VNET_LINALG_UTILS_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

namespace {
    template<typename T>
    T min_val(T a, T b) {
        return a < b ? a : b;
    }

    template<typename T>
    T max_val(T a, T b) {
        return a > b ? a : b;
    }
}

#endif //VNET_LINALG_UTILS_H
