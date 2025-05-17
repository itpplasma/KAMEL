#include "namelist.h"

#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    int a;
    double b;
    char c;
    read_namelist_unit_test(&a, &b, &c);
    assert(a == 1);
    assert(std::abs(b - 2.5) < 1e-14);
    assert(c == 'k');
    std::cout << "a: " << a << ", b: " << b << ", c: " << c << '\n';
    return 0;
}
