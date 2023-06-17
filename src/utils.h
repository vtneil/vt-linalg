#ifndef VNET_LINALG_UTILS_H
#define VNET_LINALG_UTILS_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

namespace vt {
    template<typename T>
    constexpr const T &min_val(const T &a, const T &b) { return (a < b) ? a : b; }

    template<typename T>
    constexpr const T &max_val(const T &a, const T &b) { return (a > b) ? a : b; }

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

        pair(T1 &&first, T2 &&second) : first(first), second(second) {}

        pair &operator=(const pair &other) {
            *this = pair(other);
            return *this;
        }
    };

    /**
     * Mimic std::initializer_list entirely.
     * When std::initializer_list is not available, you can
     * use this class by putting "CUSTOM_INITIALIZER_LIST" in your code.
     *
     * @tparam T
     */
    template<typename T>
    class initializer_list {
    private:
        const T *arr_;
        size_t size_;

        initializer_list(const T *array, size_t size) : arr_{array}, size_{size} {}

    public:
        initializer_list() : arr_(nullptr), size_(0) {}

        template<size_t N>
        explicit initializer_list(const T(&array)[N]) : arr_(array), size_(N) {}

        size_t size() const { return size_; }

        const T *begin() const { return arr_; }

        const T *end() const { return begin() + size(); }
    };

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
    constexpr remove_reference_t<T> &&move(T &&t) noexcept {
        return static_cast<remove_reference_t<T> &&>(t);
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

    template<typename T>
    class Iterator {
    private:
        T *ptr_;

    public:
        explicit Iterator(T *ptr) : ptr_(ptr) {}

        Iterator(const Iterator &other) : ptr_(other.ptr_) {}

        T &operator*() { return *ptr_; }

        T *operator->() { return ptr_; }

        Iterator operator+(size_t n) { return Iterator(ptr_ + n); }

        Iterator &operator++() {
            ++ptr_;
            return *this;
        }

        const Iterator operator++(int) {
            Iterator tmp(*this);
            operator++();
            return tmp;
        }

        Iterator operator-(size_t n) { return Iterator(ptr_ - n); }

        size_t operator-(Iterator other) { return ptr_ - other.ptr_; }

        Iterator &operator--() {
            --ptr_;
            return *this;
        }

        const Iterator operator--(int) {
            Iterator tmp(*this);
            operator--();
            return tmp;
        }

        bool operator==(const Iterator &other) const { return ptr_ == other.ptr_; }

        bool operator!=(const Iterator &other) const { return !operator==(other); }

        bool operator==(T *ptr) const { return ptr_ == ptr; }

        bool operator!=(T *ptr) const { return !operator==(ptr); }
    };

    template<typename T>
    bool operator==(T *ptr, const Iterator<T> &iterator) { return iterator.operator==(ptr); }

    template<typename T>
    bool operator!=(T *ptr, const Iterator<T> &iterator) { return !operator==(ptr, iterator); }
}

#endif //VNET_LINALG_UTILS_H
