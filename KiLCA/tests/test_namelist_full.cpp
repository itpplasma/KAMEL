#include "namelist.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>


bool approx(double a, double b) {
    constexpr double epsilon = 1e-14;
    return std::abs(a - b) < epsilon;
}


int main() {
    double rtor;
    double rp;
    double B0;
    char path2profiles_buffer[256];
    int calc_back;
    char flag_back;
    int N;
    double V_gal_sys;
    double V_scale;
    double m_i;
    double zele;
    double zion;
    bool flag_debug;

    read_namelist(&rtor, &rp, &B0, path2profiles_buffer, &calc_back, &flag_back, &N,
                            &V_gal_sys, &V_scale, &m_i, &zele, &zion, &flag_debug);
    std::string path2profiles{path2profiles_buffer};

    assert(approx(rtor, 165.0));
    assert(approx(rp, 67.0));
    assert(approx(B0, -17977.413));
    assert(path2profiles == "./profiles/");
    assert(calc_back == 1);
    assert(flag_back == 'f');
    assert(N == 9);
    assert(approx(V_gal_sys, -1e8));
    assert(approx(V_scale, 1.0));
    assert(approx(m_i, 2.0));
    assert(approx(zele, 1.0));
    assert(approx(zion, 1.0));
    assert(flag_debug == false);

    std::cout << "rtor: " << rtor << ", rp: " << rp << ", B0: " << B0
              << ", path2profiles: " << path2profiles
              << ", calc_back: " << calc_back
              << ", flag_back: " << flag_back
              << ", N: " << N
              << ", V_gal_sys: " << V_gal_sys
              << ", V_scale: " << V_scale
              << ", m_i: " << m_i
              << ", zele: " << zele
              << ", zion: " << zion
              << ", flag_debug: " << flag_debug
              << '\n';

    return 0;
}
