#ifndef VT_LINALG_TIE_OBJECT_H
#define VT_LINALG_TIE_OBJECT_H

#include "iterator.h"

namespace vt {
    template<typename... Ts>
    class tie_object;

    template<>
    class tie_object<> {
    };

    template<typename T, typename... Ts>
    class tie_object<T, Ts...> : private tie_object<Ts...> {
    private:
        T value;

    public:
        constexpr tie_object() = default;

        constexpr tie_object(const tie_object &) = default;

        constexpr tie_object(tie_object &&) noexcept = default;

        explicit constexpr tie_object(T value, Ts... values)
            : tie_object<Ts...>(values...), value{value} {}

        T &get() { return value; }

        tie_object<Ts...> &next() { return *this; }

        template<typename U>
        tie_object<Ts...> &operator>>(U &var) {
            var = value;
            return next();
        }

        template<typename U>
        tie_object<Ts...> &operator<<(const U &var) {
            value = var;
            return next();
        }
    };

    namespace detail {
        template<size_t, typename>
        struct tie_object_element;

        template<size_t I, typename T, typename... Ts>
        struct tie_object_element<I, tie_object<T, Ts...>> : tie_object_element<I - 1, tie_object<Ts...>> {
        };

        template<typename T, typename... Ts>
        struct tie_object_element<0, tie_object<T, Ts...>> {
            typedef T type;
        };
    }  // namespace detail

    template<size_t I, typename T, typename... Ts>
    vt::enable_if_t<(I == 0), T> get_tie(tie_object<T, Ts...> &t) {
        return t.get();
    }

    template<size_t I, typename T, typename... Ts>
    vt::enable_if_t<(I != 0), typename detail::tie_object_element<I, tie_object<T, Ts...>>::type> get_tie(tie_object<T, Ts...> &t) {
        return get_tie<I - 1>(t.next());
    }

    template<typename... Ts>
    tie_object<Ts...> make_tie_object(Ts... ts) {
        return tie_object<Ts...>(ts...);
    }

    template<typename... Ts>
    constexpr tie_object<Ts &...> tie(Ts &...ts) noexcept {
        return tie_object<Ts &...>(ts...);
    }

    template<typename... Ts>
    constexpr tie_object<const Ts &...> tie(const Ts &...ts) noexcept {
        return tie_object<const Ts &...>(ts...);
    }
}  // namespace vt

#endif  //VT_LINALG_TIE_OBJECT_H
