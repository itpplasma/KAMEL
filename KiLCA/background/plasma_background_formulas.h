#ifndef KILCA_PLASMA_BACKGROUND_FORMULAS_H
#define KILCA_PLASMA_BACKGROUND_FORMULAS_H

struct kilca_background_values {
    double lee;
    double lei;
    double lii;
    double nuee;
    double nuei;
    double nuie;
    double nuii;
    double nue;
    double nui;
    double vti;
    double vte;
};

kilca_background_values evaluate_kilca_background(
    double density, double ion_temperature, double electron_temperature,
    double ion_charge_number, double ion_mass_number, double ion_mass,
    double electron_mass, double velocity_factor);

#endif
