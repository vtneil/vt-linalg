#ifndef VNET_LINALG_ITERATOR_H
#define VNET_LINALG_ITERATOR_H

#include "standard_utility.h"

namespace vt {
    template<typename T>
    class iterator {
    private:
        T *ptr_;

    public:
        explicit iterator(T *ptr) : ptr_(ptr) {}

        iterator(const iterator &other) : ptr_(other.ptr_) {}

        T &operator*() { return *ptr_; }

        T *operator->() { return ptr_; }

        iterator operator+(size_t n) { return iterator(ptr_ + n); }

        iterator &operator++() {
            ++ptr_;
            return *this;
        }

        const iterator operator++(int) {
            iterator tmp(*this);
            operator++();
            return tmp;
        }

        iterator operator-(size_t n) { return iterator(ptr_ - n); }

        size_t operator-(iterator other) { return ptr_ - other.ptr_; }

        iterator &operator--() {
            --ptr_;
            return *this;
        }

        const iterator operator--(int) {
            iterator tmp(*this);
            operator--();
            return tmp;
        }

        bool operator==(const iterator &other) const { return ptr_ == other.ptr_; }

        bool operator!=(const iterator &other) const { return !operator==(other); }

        bool operator==(T *ptr) const { return ptr_ == ptr; }

        bool operator!=(T *ptr) const { return !operator==(ptr); }
    };

    template<typename T>
    bool operator==(T *ptr, const iterator<T> &iterator) { return iterator.operator==(ptr); }

    template<typename T>
    bool operator!=(T *ptr, const iterator<T> &iterator) { return !operator==(ptr, iterator); }
}

#endif //VNET_LINALG_ITERATOR_H
