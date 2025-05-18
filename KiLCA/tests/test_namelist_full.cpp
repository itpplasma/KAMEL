#include "test_namelist.hpp"

#include <iostream>
#include <string>
#include <vector>


int main() {
    // background
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

    // eigmode
    char fname_buffer[256];
    bool search_flag;
    int rdim;
    double rfmin;
    double rfmax;
    int idim;
    double ifmin;
    double ifmax;
    bool stop_flag;
    double eps_res;
    double eps_abs;
    double eps_rel;
    double delta;
    bool test_roots;
    int Nguess;
    int kmin;
    int kmax;
    std::vector<std::complex<double>> fstart(100);

    // debuggroup
    bool flag_debug;

    read_namelist(&rtor, &rp, &B0, path2profiles_buffer, &calc_back, &flag_back, &N, &V_gal_sys,
                  &V_scale, &m_i, &zele, &zion,
                  fname_buffer, &search_flag, &rdim, &rfmin, &rfmax, &idim, &ifmin, &ifmax,
                  &stop_flag, &eps_res, &eps_abs, &eps_rel, &delta, &test_roots, &Nguess,
                  &kmin, &kmax, fstart.data(),
                  &flag_debug);

    // convert buffers
    std::string path2profiles{path2profiles_buffer};
    std::string fname{fname_buffer};
    fstart.resize(Nguess);

    std::cout << "rtor: " << rtor
              << ", rp: " << rp
              << ", B0: " << B0
              << ", path2profiles: " << path2profiles
              << ", calc_back: " << calc_back
              << ", flag_back: " << flag_back
              << ", N: " << N
              << ", V_gal_sys: " << V_gal_sys
              << ", V_scale: " << V_scale
              << ", m_i: " << m_i
              << ", zele: " << zele
              << ", zion: " << zion
              << ", fname: " << fname
              << ", search_flag: " << search_flag
              << ", rdim: " << rdim
              << ", rfmin: " << rfmin
              << ", rfmax: " << rfmax
              << ", idim: " << idim
              << ", ifmin: " << ifmin
              << ", ifmax: " << ifmax
              << ", stop_flag: " << stop_flag
              << ", eps_res: " << eps_res
              << ", eps_abs: " << eps_abs
              << ", eps_rel: " << eps_rel
              << ", delta: " << delta
              << ", test_roots: " << test_roots
              << ", Nguess: " << Nguess
              << ", kmin: " << kmin
              << ", kmax: " << kmax
              << ", fstart: [";
    for (size_t i = 0; i < fstart.size(); ++i) {
        std::cout << fstart[i];
        if (i < fstart.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]"
              << ", flag_debug: " << flag_debug
              << '\n';

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

    return 0;
}
