#ifndef VT_LINALG_INITIALIZER_LIST_H
#define VT_LINALG_INITIALIZER_LIST_H

#include "standard_utility.h"

namespace vt {
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
}

#endif //VT_LINALG_INITIALIZER_LIST_H
