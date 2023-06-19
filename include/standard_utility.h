#ifndef VNET_LINALG_STANDARD_UTILITY_H
#define VNET_LINALG_STANDARD_UTILITY_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define __VT_FORCE_INLINE __attribute__((always_inline))

namespace vt {
    typedef double real_t;

    template<typename T>
    __VT_FORCE_INLINE constexpr const T &min(const T &a, const T &b) { return (a < b) ? a : b; }

    template<typename T, typename... Ts>
    __VT_FORCE_INLINE constexpr const T &min(const T &a, const T &b, const Ts &... args) {
        return min(min(a, b), args...);
    }

    template<typename T>
    __VT_FORCE_INLINE constexpr const T &max(const T &a, const T &b) { return (a > b) ? a : b; }

    template<typename T, typename... Ts>
    __VT_FORCE_INLINE constexpr const T &max(const T &a, const T &b, const Ts &... args) {
        return max(max(a, b), args...);
    }

    template<typename T>
    struct remove_reference {
        using type = T;
    };

    template<typename T>
    struct remove_reference<T &> {
        using type = T;
    };

    template<typename T>
    struct remove_reference<T &&> {
        using type = T;
    };

    template<typename T>
    using remove_reference_t = typename vt::remove_reference<T>::type;

    /**
     * Mimic std::move.
     *
     * @tparam T
     * @param t
     * @return
     */
    template<typename T>
    __VT_FORCE_INLINE constexpr
    remove_reference_t<T> &&move(T &&t) noexcept { return static_cast<remove_reference_t<T> &&>(t); }

    /**
     * Mimic std::swap
     *
     * @tparam T
     * @param a
     * @param b
     */
    template<typename T>
    void swap(T &a, T &b) {
        T tmp = vt::move(a);
        a = vt::move(b);
        b = vt::move(tmp);
    }

    /**
     * Mimic std::is_same.
     *
     * @tparam T
     * @tparam U
     */
    template<typename T, typename U>
    struct is_same {
        static constexpr bool value = false;
    };

    /**
     * Mimic std::is_same.
     *
     * @tparam T
     */
    template<typename T>
    struct is_same<T, T> {
        static constexpr bool value = true;
    };

    /**
     * Mimic std::fill
     *
     * @tparam ForwardIt
     * @tparam T
     * @param first
     * @param last
     * @param value
     */
    template<typename ForwardIt, typename T>
    void fill(ForwardIt first, ForwardIt last, const T &value) {
        for (; first != last; static_cast<void>(++first)) *first = value;
    }

    /**
     * Mimic std::copy
     *
     * @tparam InputIt
     * @tparam OutputIt
     * @param first
     * @param last
     * @param d_first
     * @return
     */
    template<typename InputIt, typename OutputIt>
    OutputIt copy(InputIt first, InputIt last, OutputIt d_first) {
        for (; first != last; static_cast<void>(++first), static_cast<void>(++d_first)) *d_first = *first;
        return d_first;
    }
}

#endif //VNET_LINALG_STANDARD_UTILITY_H
