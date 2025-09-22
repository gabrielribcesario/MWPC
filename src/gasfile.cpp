#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <filesystem>

#include <Garfield/MediumMagboltz.hh>
#include <Garfield/FundamentalConstants.hh>

using namespace Garfield;

int main(int argc, char **argv) {
    // Number of collisions (in multiples of 10^7) over which the electron is traced by Magboltz
    // Recommended value is nCollisions >= 10
    const int nCollisions = 10;
    // Number of electric field points to be evaluated
    const size_t nE = 40;
    // Sets the electric field range [V/cm] to be covered by the gas table
    const double Emin = 100, Emax = 1.E6;
    // Use logarithmic spacing
    const bool useLog = true;

    // Gas pressure [Torr].
    const double gas_pressure = 3. * AtmosphericPressure;
    // Gas temperature [K]
    const double gas_temperature = 293.15;

    // Saves the gas files here
    const std::string gas_folder = "../data/gas/1atm";
    std::filesystem::create_directories(gas_folder);

    // Sets the gas composition, temperature and pressure
    MediumMagboltz gas("Ar", 80., "CO2", 20.);
    gas.SetTemperature(gas_temperature);
    gas.SetPressure(gas_pressure);
    // Sets the points of evaluation
    gas.SetFieldGrid(Emin, Emax, nE, useLog); 
    // Runs the Magboltz algorithm and generates the gas table
    gas.GenerateGasTable(nCollisions, true);

    // Parse .gas file name
    std::vector<std::string> gases(6, "");
    std::vector<double> fracs(6, 0.0);
    gas.GetComposition(gases[0], fracs[0], 
                       gases[1], fracs[1], 
                       gases[2], fracs[2],
                       gases[3], fracs[3], 
                       gases[4], fracs[4], 
                       gases[5], fracs[5]);
    std::string gas_name = gases[0] + "_" + std::to_string((unsigned int) std::round(fracs[0] * 100.0));
    for (unsigned int i = 1; !gases[i].empty() && i < 6; ++i) {
        gas_name += "_" + gases[i] + "_" + std::to_string((unsigned int) std::round(fracs[i] * 100.0));
    }

    // If a gas table already exists, merge both previous and current files
    std::string gasfile = gas_folder + "/" + gas_name + ".gas";
    if(std::filesystem::exists(gasfile)) {
        if (!gas.MergeGasFile(gasfile, true)) { gasfile = gas_folder + "/new_" + gas_name + ".gas"; }
    }
    // Saves the .gas file
    gas.WriteGasFile(gasfile);

    return EXIT_SUCCESS;
}