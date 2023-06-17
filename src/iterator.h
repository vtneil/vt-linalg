#ifndef VNET_LINALG_ITERATOR_H
#define VNET_LINALG_ITERATOR_H

#include "standard_utility.h"

namespace vt {
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

#endif //VNET_LINALG_ITERATOR_H
