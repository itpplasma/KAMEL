#include "calc_flre_quants.h"

#include <cmath>
#include <complex>
#include <cstdio>

namespace {
bool close(double actual, double expected, double tolerance = 1.0e-14) {
    return std::abs(actual - expected) <= tolerance;
}

int check(const char* name, double actual, double expected, double tolerance = 1.0e-14) {
    if (close(actual, expected, tolerance)) {
        return 0;
    }

    std::fprintf(stderr, "%s: expected %.17e, got %.17e (difference %.17e)\n", name, expected,
        actual, actual - expected);
    return 1;
}
} // namespace

int main() {
    int failures = 0;

    const std::complex<double> density(1.0, 2.0);
    const std::complex<double> electric_field[3] = {std::complex<double>(3.0, 4.0),
        std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0)};
    const std::complex<double> current_density[3] = {std::complex<double>(1.0, 1.0),
        std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0)};
    const std::complex<double> magnetic_field[3] = {std::complex<double>(0.0, 0.0),
        std::complex<double>(0.0, 0.0), std::complex<double>(2.0, 3.0)};
    double force[3];

    calc_time_averaged_lorentz_force(
        2.0, density, electric_field, current_density, magnetic_field, force);

    failures += check("electrostatic force", force[0], 11.0);
    failures += check("magnetic cross-product sign", force[1], -2.5 / c);
    failures += check("zero force component", force[2], 0.0);

    double torque[3];
    calc_cylindrical_torque_density(5.0, 7.0, force, torque);

    failures += check("poloidal torque moment arm", torque[1], -12.5 / c);
    failures += check("toroidal torque moment arm", torque[2], 0.0);

    double radius[3] = {1.0, 2.0, 4.0};
    double density_profile[3] = {1.0, -2.0, 3.0};
    double integrated[3];

    integrate_over_cylinder(3, radius, density_profile, 10.0, integrated);

    failures += check("integral origin", integrated[0], 0.0);
    failures += check("integral first interval", integrated[1], -15.0);
    failures += check("integral signed total", integrated[2], 65.0);

    if (failures != 0) {
        return 1;
    }

    std::puts("KiLCA Lorentz torque identities passed");
    return 0;
}
