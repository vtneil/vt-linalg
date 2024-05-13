#ifndef VT_LINALG_ITERATOR_H
#define VT_LINALG_ITERATOR_H

#include "standard_utility.h"

namespace vt {
    template<typename T>
    class iterator {
    private:
        T *ptr_;

    public:
        explicit constexpr iterator(T *ptr) : ptr_(ptr) {}

        constexpr iterator(const iterator &other) : ptr_(other.ptr_) {}

        T &operator*() { return *ptr_; }

        T *operator->() { return ptr_; }

        constexpr iterator operator+(size_t n) const { return iterator(ptr_ + n); }

        iterator &operator++() {
            ++ptr_;
            return *this;
        }

        const iterator operator++(int) {
            iterator tmp(*this);
            operator++();
            return tmp;
        }

        constexpr iterator operator-(size_t n) const { return iterator(ptr_ - n); }

        constexpr size_t operator-(iterator other) const { return ptr_ - other.ptr_; }

        iterator &operator--() {
            --ptr_;
            return *this;
        }

        const iterator operator--(int) {
            iterator tmp(*this);
            operator--();
            return tmp;
        }

        constexpr bool operator==(const iterator &other) const { return ptr_ == other.ptr_; }

        constexpr bool operator!=(const iterator &other) const { return !operator==(other); }

        constexpr bool operator==(T *ptr) const { return ptr_ == ptr; }

        constexpr bool operator!=(T *ptr) const { return !operator==(ptr); }
    };

    template<typename T>
    constexpr bool operator==(T *ptr, const iterator<T> &iterator) { return iterator.operator==(ptr); }

    template<typename T>
    constexpr bool operator!=(T *ptr, const iterator<T> &iterator) { return !operator==(ptr, iterator); }

    template<typename T>
    using const_iterator = iterator<const T>;

    namespace detail {
        template<typename T, size_t Size, size_t Index>
        struct output_assignment_chain {
            const T (&arr_ref)[Size];

            explicit output_assignment_chain(const T (&arr)[Size]) : arr_ref{arr} {}

            static_assert(Size > 0, "Size must be greater than 0.");
            using next_iterator = output_assignment_chain<T, Size, Index + 1>;

            FORCE_INLINE next_iterator operator>>(T &var) const {
                var = arr_ref[Index];
                static_assert(Index < Size, "Can\'t assign index greater than or equal to size.");
                return next_iterator(arr_ref);
            }
        };

        template<typename T, size_t Size, size_t Index>
        struct input_assignment_chain {
            T (&arr_ref)
            [Size];

            explicit input_assignment_chain(T (&arr)[Size]) : arr_ref{arr} {}

            static_assert(Size > 0, "Size must be greater than 0.");
            using next_iterator = input_assignment_chain<T, Size, Index + 1>;

            FORCE_INLINE next_iterator operator<<(const T &var) {
                arr_ref[Index] = var;
                static_assert(Index < Size, "Can\'t assign index greater than or equal to size.");
                return next_iterator(arr_ref);
            }
        };
    }// namespace detail
}// namespace vt

#endif//VT_LINALG_ITERATOR_H
