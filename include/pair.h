#ifndef VT_LINALG_PAIR_H
#define VT_LINALG_PAIR_H

#include "standard_utility.h"

namespace vt {
    /**
     * Mimic std::pair.
     *
     * @tparam T1
     * @tparam T2
     */
    template<typename T1, typename T2>
    class pair {
    public:
        T1 first;
        T2 second;

        pair() = default;

        pair(const pair &other) : first(other.first), second(other.second) {}

        pair(const T1 &first, const T2 &second) : first(first), second(second) {}

        pair(T1 &&first, T2 &&second) noexcept : first(move(first)), second(move(second)) {}

        pair &operator=(const pair &other) {
            *this = pair(other);
            return *this;
        }
    };

    template<typename T1, typename T2>
    constexpr pair<T1, T2> make_pair(const T1 &t1, const T2 &t2) { return pair<T1, T2>(t1, t2); }

    template<typename T1, typename T2>
    constexpr pair<T1, T2> make_pair(T1 &&t1, T2 &&t2) { return pair<T1, T2>(t1, t2); }
}// namespace vt

#endif//VT_LINALG_PAIR_H
