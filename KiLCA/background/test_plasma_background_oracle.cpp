#include "plasma_background_formulas.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace {

bool close_value(const char* label, double actual, double expected)
{
    const double error = std::abs(actual - expected) / std::max(1.0, std::abs(expected));
    if (error <= 5.0e-12) return true;
    std::cerr << "FAIL: " << label << " actual=" << actual
              << " expected=" << expected << " scaled_error=" << error << '\n';
    return false;
}

} // namespace

int main()
{
    std::ifstream oracle("plasma_background.dat");
    if (!oracle) {
        std::cerr << "FAIL: cannot open plasma_background.dat\n";
        return 1;
    }

    std::string line;
    int matched_rows = 0;
    bool passed = true;
    while (std::getline(oracle, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream input(line);
        int kind = 0, case_id = 0, species = 0, reserved = 0;
        double value[20] = {};
        input >> kind >> case_id >> species >> reserved;
        for (double& item : value) input >> item;
        if (!input) {
            std::cerr << "FAIL: malformed plasma-background oracle row\n";
            return 1;
        }
        if (kind != 4) continue;

        const kilca_background_values actual = evaluate_kilca_background(
            value[0], value[1], value[2], value[3], value[4], value[5],
            value[6], value[7]);
        const double result[] = {actual.lee, actual.lei, actual.lii,
                                 actual.nuee, actual.nuei, actual.nuie,
                                 actual.nuii, actual.nue, actual.nui,
                                 actual.vti, actual.vte};
        const char* label[] = {"Lee", "Lei", "Lii", "nuee", "nuei", "nuie",
                               "nuii", "nue", "nui", "Vti", "Vte"};
        for (int i = 0; i < 11; ++i) {
            passed = close_value(label[i], result[i], value[8 + i]) && passed;
        }
        ++matched_rows;
    }

    if (matched_rows != 1) {
        std::cerr << "FAIL: expected one KiLCA background oracle row, got "
                  << matched_rows << '\n';
        return 1;
    }
    if (!passed) return 1;
    std::cout << "PASS: KiLCA legacy plasma-background oracle\n";
    return 0;
}
