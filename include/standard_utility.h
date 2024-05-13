#ifndef VT_LINALG_STANDARD_UTILITY_H
#define VT_LINALG_STANDARD_UTILITY_H

#include "standard_constants.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef FORCE_INLINE
#undef FORCE_INLINE
#endif
#define FORCE_INLINE __attribute__((always_inline))

#ifdef NO_INLINE
#undef NO_INLINE
#endif
#define NO_INLINE __attribute__((noinline))

namespace vt {
    using real_t = double;

    template<typename T>
    FORCE_INLINE constexpr const T &min(const T &a, const T &b) { return (a < b) ? a : b; }

    template<typename T, typename... Ts>
    FORCE_INLINE constexpr const T &min(const T &a, const T &b, const Ts &...args) {
        return min(min(a, b), args...);
    }

    template<typename T>
    FORCE_INLINE constexpr const T &max(const T &a, const T &b) { return (a > b) ? a : b; }

    template<typename T, typename... Ts>
    FORCE_INLINE constexpr const T &max(const T &a, const T &b, const Ts &...args) {
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
    FORCE_INLINE constexpr remove_reference_t<T> &&move(T &&t) noexcept { return static_cast<remove_reference_t<T> &&>(t); }

    template<class T>
    struct remove_const {
        typedef T type;
    };
    template<class T>
    struct remove_const<const T> {
        typedef T type;
    };

    template<typename T>
    using remove_const_t = typename vt::remove_const<T>::type;

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
        a     = vt::move(b);
        b     = vt::move(tmp);
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

    namespace detail {
        template<typename>
        struct always_false {
            static constexpr bool value = false;
        };

        template<typename T>
        struct type_identity {
            using type = T;
        };

        template<typename T>
        auto try_add_lvalue_reference(int) -> type_identity<T &>;

        template<typename T>
        auto try_add_lvalue_reference(...) -> type_identity<T>;

        template<typename T>
        auto try_add_rvalue_reference(int) -> type_identity<T &&>;

        template<typename T>
        auto try_add_rvalue_reference(...) -> type_identity<T>;
    }  // namespace detail

    template<typename T>
    struct add_lvalue_reference : decltype(vt::detail::try_add_lvalue_reference<T>(0)) {
    };

    template<typename T>
    struct add_rvalue_reference : decltype(vt::detail::try_add_rvalue_reference<T>(0)) {
    };

    template<typename T>
    typename vt::add_rvalue_reference<T>::type declval() noexcept {
        static_assert(vt::detail::always_false<T>::value, "declval is not allowed in an evaluated context.");
    }

    template<typename T, T v>
    struct integral_constant {
        static constexpr const T value = v;
        typedef T value_type;
        typedef integral_constant type;

        explicit constexpr operator value_type() const noexcept { return value; }
    };

    typedef vt::integral_constant<bool, true> true_type;

    typedef vt::integral_constant<bool, false> false_type;

    template<typename T>
    struct is_lvalue_reference : public vt::false_type {
    };

    template<typename T>
    struct is_lvalue_reference<T &> : public vt::true_type {
    };

    template<typename T>
    struct is_rvalue_reference : public vt::false_type {
    };

    template<typename T>
    struct is_rvalue_reference<T &&> : public vt::true_type {
    };

    template<typename T>
    struct is_reference : public vt::false_type {
    };

    template<typename T>
    struct is_reference<T &> : public vt::true_type {
    };

    template<typename T>
    struct is_reference<T &&> : public vt::true_type {
    };

    template<typename T>
    inline constexpr T &&forward(vt::remove_reference_t<T> &t) noexcept {
        return static_cast<T &&>(t);
    }

    template<typename T>
    inline constexpr T &&forward(vt::remove_reference_t<T> &&t) noexcept {
        static_assert(!is_lvalue_reference<T>::value, "cannot forward an rvalue as an lvalue");
        return static_cast<T &&>(t);
    }

    template<bool ctime_condition>
    struct if_constexpr;

    template<>
    struct if_constexpr<true> {
        template<typename Func, typename... Args>
        static void run(Func &&func, Args &&...args) { func(vt::forward<Args>(args)...); }
    };

    template<>
    struct if_constexpr<false> {
        template<typename Func, typename... Args>
        static void run(Func &&, Args &&...) {}
    };

    template<size_t... Is>
    struct index_sequence {
    };

    template<size_t N, size_t... Is>
    struct make_index_sequence : public make_index_sequence<N - 1, N - 1, Is...> {
    };

    template<size_t... Is>
    struct make_index_sequence<0, Is...> : public index_sequence<Is...> {
    };

    template<bool B, typename T = void>
    struct enable_if {
    };

    template<typename T>
    struct enable_if<true, T> {
        typedef T type;
    };

    template<bool B, typename T = void>
    using enable_if_t = typename enable_if<B, T>::type;

    template<bool condition, typename T>
    constexpr T choose_if(const T &value_if_true, const T &value_if_false) {
        return condition ? value_if_true : value_if_false;
    }

    namespace detail {
        template<typename T, size_t N>
        struct integral_coefficient_helper {
            static constexpr real_t value =
                    (1.0 / static_cast<real_t>(N + 1)) * integral_coefficient_helper<T, N - 1>::value;
        };

        template<typename T>
        struct integral_coefficient_helper<T, 1> {
            static constexpr real_t value = 0.5;
        };
    }  // namespace detail

    template<size_t N>
    constexpr real_t integral_coefficient() {
        return detail::integral_coefficient_helper<real_t, N>::value;
    }
}  // namespace vt

#endif  //VT_LINALG_STANDARD_UTILITY_H
