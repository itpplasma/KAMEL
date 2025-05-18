#include "test_namelist.hpp"

#include <iostream>
#include <string>

int main() {
    int a;
    double b;
    char c;
    char d_buffer[16];
    int e[3];
    std::complex<double> f[3];

    read_namelist_unit_test(&a, &b, &c, d_buffer, e, f);
    std::string d{d_buffer};

    std::cout << "a: " << a << ", b: " << b << ", c: " << c << ", d: " << d
              << ", e: [" << e[0] << ", " << e[1] << ", " << e[2] << "]"
              << ", f: [" << f[0] << ", " << f[1] << ", " << f[2] << "]"
              << '\n';

    assert(a == 1);
    assert(approx(b, 2.5));
    assert(c == 'k');
    assert(d == "hello");
    assert(e[0] == 1);
    assert(e[1] == 2);
    assert(e[2] == 3);

    return 0;
}
