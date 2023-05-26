#ifndef VNET_LINALG_UTILS_H
#define VNET_LINALG_UTILS_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


template<typename T>
inline T min_val(T a, T b) { return a < b ? a : b; }

template<typename T>
inline T max_val(T a, T b) { return a > b ? a : b; }


template<typename T1, typename T2>
class Pair {
public:
    T1 first;
    T2 second;

    Pair() = default;

    Pair(const Pair &other) : first(other.first), second(other.second) {}

    Pair(const T1 &first, const T2 &second) : first(first), second(second) {}

    Pair(T1 &&first, T2 &&second) : first(first), second(second) {}

    Pair &operator=(const Pair &other) {
        *this = Pair(other);
        return *this;
    }

    T1 &lower() { return first; }

    T1 &upper() { return second; }
};

#endif //VNET_LINALG_UTILS_H
