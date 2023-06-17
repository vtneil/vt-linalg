#include <iostream>
#include "../src/standard_utility.h"
#include "../src/buffer.h"

using namespace vt;

template<typename T, size_t N>
void print_buffer(buffer_t<T, N> buffer) {
    std::cout << "Capacity = " << buffer.capacity() << '\n';
    std::cout << "Size = " << buffer.size() << '\n';
    std::cout << "Start Index = " << (&buffer.front() - buffer.arr()) << " has content = " << buffer.front() << '\n';
    std::cout << "Content = ";
    for (size_t i = 0; i < buffer.capacity(); ++i) std::cout << buffer.arr()[i] << ' ';
    std::cout << "\n\n";
}

int main() {
    std::cout << vt::max(1, 2, 3, 4, 5, 6) << '\n';
    buffer_t<int, 8> buffer = make_buffer<int, 8>();
    int x = 1;

    print_buffer(buffer);

    buffer.push(x);
    buffer.push(5);

    print_buffer(buffer);

    buffer.pop();
    buffer.pop();

    print_buffer(buffer);

    buffer.push(8);
    buffer.push(9);
    buffer.push(10);
    buffer.push(11);
    buffer.push(12);
    buffer.push(13);
    buffer.push(14);

    print_buffer(buffer);

    while (!buffer.empty()) {
        std::cout << buffer.front() << ' ';
        buffer.pop();
    }

    return 0;
}
