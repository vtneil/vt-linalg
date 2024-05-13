#ifndef VT_LINALG_CIRCULAR_BUFFER_H
#define VT_LINALG_CIRCULAR_BUFFER_H

#include "standard_utility.h"

namespace vt {
    namespace {
        template<size_t Capacity>
        constexpr size_t cyclic(size_t current) { return current % Capacity; }
    }  // namespace

    template<typename T, size_t Capacity>
    class circular_buffer_static_t {
    public:
        static_assert(Capacity > 0, "Buffer capacity must not be 0.");

    private:
        T arr_[Capacity]    = {};
        size_t start_index_ = 0;
        size_t size_        = 0;

    public:
        circular_buffer_static_t() = default;

        circular_buffer_static_t(const circular_buffer_static_t &other) {
            start_index_ = other.start_index_;
            size_        = other.size_;
            vt::copy(other.arr_, other.arr_ + Capacity, arr_);
        }

        circular_buffer_static_t(circular_buffer_static_t &&other) noexcept {
            start_index_       = other.start_index_;
            size_              = other.size_;
            other.start_index_ = 0;
            other.size_        = 0;
            for (size_t i = 0; i < Capacity; ++i) arr_[i] = vt::move(other.arr_[i]);
        }

        void push(const T &t) {
            if (full()) return;
            arr_[cyclic<Capacity>(start_index_ + size_)] = move(t);
            static_cast<void>(++size_);
        }

        void push(T &&t) {
            if (full()) return;
            arr_[cyclic<Capacity>(start_index_ + size_)] = move(t);
            static_cast<void>(++size_);
        }

        void pop() {
            if (empty()) return;
            start_index_ = cyclic<Capacity>(start_index_ + 1);
            static_cast<void>(--size_);
        }

        T &front() { return arr_[start_index_]; }

        constexpr T &front() const { return arr_[start_index_]; }

        T &back() { return arr_[cyclic<Capacity>(start_index_ + size_ - 1)]; }

        constexpr T &back() const { return arr_[cyclic<Capacity>(start_index_ + size_ - 1)]; }

        constexpr bool empty() const { return (size_ == 0); }

        constexpr bool full() const { return !available_for(1); }

        constexpr bool available_for(size_t size) const { return (Capacity - size_ >= size); }

        constexpr size_t size() const { return size_; }

        constexpr size_t capacity() const { return Capacity; }

        constexpr const T *arr() const { return arr_; }

        constexpr circular_buffer_static_t copy() const { return circular_buffer_static_t(*this); }
    };

    class circular_buffer_dynamic_t {
    public:
        circular_buffer_dynamic_t() = delete;
    };

    template<typename T, size_t Capacity>
    using buffer_t = circular_buffer_static_t<T, Capacity>;

    template<typename T, size_t Capacity>
    constexpr buffer_t<T, Capacity> make_buffer() { return buffer_t<T, Capacity>(); }
}  // namespace vt

#endif  //VT_LINALG_CIRCULAR_BUFFER_H
