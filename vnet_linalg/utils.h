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

template<typename T>
inline void swap_val(T &a, T &b) {
    T tmp = a;
    a = b;
    b = tmp;
}

template<typename T>
inline void swap_val(T *&a, T *&b) {
    T *tmp = a;
    a = b;
    b = tmp;
}


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

template<typename T>
class InitializerList {
private:
    const T *arr_;
    size_t size_;

    InitializerList(const T *array, size_t size) : arr_{array}, size_{size} {}

public:
    InitializerList() : arr_(nullptr), size_(0) {}

    template<size_t N>
    explicit InitializerList(const T(&array)[N]) : arr_(array), size_(N) {}

    size_t size() const { return size_; }

    const T *begin() const { return arr_; }

    const T *end() const { return begin() + size(); }
};

#endif //VNET_LINALG_UTILS_H
