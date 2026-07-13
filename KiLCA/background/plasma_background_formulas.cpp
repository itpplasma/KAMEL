#include "plasma_background_formulas.h"

#include "constants.h"

#include <cmath>

kilca_background_values evaluate_kilca_background(
    double density, double ion_temperature, double electron_temperature,
    double ion_charge_number, double ion_mass_number, double ion_mass,
    double electron_mass, double velocity_factor)
{
    kilca_background_values value;
    value.lee = 23.5 - std::log(std::sqrt(density) /
                               std::pow(electron_temperature, 1.25)) -
                (std::log(electron_temperature) - 2.0) / 4.0;
    value.lei = 30.0 - std::log(std::sqrt(density) /
                               std::pow(ion_temperature, 1.5) *
                               std::pow(ion_charge_number, 2.0) /
                               ion_mass_number);
    value.lii = 23.0 - std::log(std::pow(ion_charge_number, 3.0) /
                               std::pow(ion_temperature, 1.5) *
                               std::sqrt(2.0 * density));

    value.nuee = 5.8e-6 * density * value.lee /
                 (std::sqrt(electron_temperature) * velocity_factor *
                  electron_temperature);
    value.nuei = 7.7e-6 * density * value.lei *
                 std::pow(ion_charge_number, 2.0) /
                 std::pow(velocity_factor * electron_temperature, 1.5);
    value.nuie = 3.2e-9 * density * value.lei *
                 std::pow(ion_charge_number, 2.0) /
                 (ion_mass_number * std::sqrt(electron_temperature) *
                  velocity_factor * ion_temperature);
    value.nuii = 1.4e-7 * density * value.lii *
                 std::pow(ion_charge_number, 4.0) /
                 (std::sqrt(ion_mass_number) * std::sqrt(ion_temperature) *
                  velocity_factor * ion_temperature);

    value.nue = value.nuee + value.nuei + 10.0;
    value.nui = value.nuie + value.nuii + 10.0;
    value.vti = std::sqrt(boltz * ion_temperature / ion_mass);
    value.vte = std::sqrt(boltz * electron_temperature / electron_mass);
    return value;
}
