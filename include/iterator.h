#ifndef VNET_LINALG_ITERATOR_H
#define VNET_LINALG_ITERATOR_H

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
}

#endif //VNET_LINALG_ITERATOR_H
