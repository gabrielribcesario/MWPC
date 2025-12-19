#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <omp.h>

#include <TFile.h>
#include <TNtupleD.h>
#include <TH1D.h>

#include <Garfield/MediumConductor.hh>
#include <Garfield/MediumMagboltz.hh>
#include <Garfield/GeometrySimple.hh>
#include <Garfield/SolidWire.hh>
#include <Garfield/ComponentNeBem3d.hh>
#include <Garfield/Sensor.hh>
#include <Garfield/AvalancheMicroscopic.hh>

#include "WirePlane.hpp"

using namespace Garfield;

int main(int argc, char **argv) {
    const char *outputFileName = "fluctuations.root";

    // Number of events to simulate
    const size_t nEvents = 1000;

    // Cathode-Anode planes spacing [cm]
    constexpr double planeGap = 0.5;

    // Wire orientation
    const char cathodeOrientation[] = "x";
    // Number of cathode wire gaps
    constexpr size_t cathodeNGaps = 8;
    // Cathode wire spacing [cm]
    constexpr double cathodeGap = 0.2;
    // Cathode wire radius [cm]
    constexpr double cathodeRadius = 0.0025;
    // Cathode wire potential [V]
    constexpr double cathodePotential = 0.;

    // Wire orientation
    const char anodeOrientation[] = "y";
    // Number of anode wire gaps
    constexpr size_t anodeNGaps = 4;
    // Anode wire spacing [cm]
    constexpr double anodeGap = 0.4;
    // Anode wire radius [cm]
    constexpr double anodeRadius = 0.001;
    // Anode wire potential [V]
    constexpr double anodePotential = 2100.;

    // Cathode wire length [cm]
    constexpr double cathodeLength = anodeNGaps * anodeGap;
    // Anode wire length [cm]
    constexpr double anodeLength = cathodeNGaps * cathodeGap;

    // Wire material
    MediumConductor metal;

    // Medium
    MediumMagboltz gas("Ar", 80., "CO2", 20.);
    gas.EnablePenningTransfer();
    gas.SetMaxElectronEnergy(200.);

    // Detector geometry
    GeometrySimple geo;
    geo.SetMedium(&gas);

    WirePlane cathodeUpperPlane(cathodeNGaps, planeGap, cathodeGap, 
                                cathodeRadius, cathodePotential, cathodeOrientation, 
                                &geo, &metal);

    WirePlane anodePlane(anodeNGaps, 0.0, anodeGap, 
                         anodeRadius, anodePotential, anodeOrientation, 
                         &geo, &metal);

    WirePlane cathodeLowerPlane(cathodeNGaps, -planeGap, cathodeGap, 
                                cathodeRadius, cathodePotential, cathodeOrientation, 
                                &geo, &metal);

    // Electric field calculation
    ComponentNeBem3d cmp;
    cmp.SetGeometry(&geo);
    cmp.SetTargetElementSize(0.1);
    cmp.UseLUInversion();
    cmp.SetNumberOfThreads(36);
    cmp.Initialise();

    // Bounding box
    constexpr double xlim = cathodeGap * cathodeNGaps * 0.5;
    constexpr double ylim = anodeGap * anodeNGaps * 0.5;
    constexpr double zlim = planeGap * 1.1;

    // Interface between the transport classes and the component
    Sensor sensor;
    sensor.AddComponent(&cmp);
    sensor.SetArea(-xlim, -ylim, -zlim,
                    xlim,  ylim,  zlim);

    // Output .root file
    TFile rootOut(outputFileName, "recreate", "", 505);

    TH1::StatOverflows(true);

    TH1D hElectrons("hElectrons", "Electron Avalanche Size", 200, 0, 2000);
    TH1D hIons("hIons", "Ion Avalanche Size", 200, 0, 2000);
    TH1D hEnergy("hEnergy", "Electron Energy Distribution", 500, 0., 50.);

    // Charge transport
    AvalancheMicroscopic aval(&sensor);
    aval.EnableRKNSteps();
    aval.SetRKNTolerance(1.e-8, 1.e-4);
    // aval.EnableElectronEnergyHistogramming(&hEnergy);

    unsigned int nIon = 0;
    unsigned int nExc = 0;

    // #pragma omp parallel for firstprivate(aval)
    for (int i = 0; i < 1; ++i) {
        // gas.ResetCollisionCounters();

        double x0 = 0.0, y0 = 0.0, z0 = planeGap * 0.9, t0 = 0.0, e0 = 1.0;
        aval.AvalancheElectron(x0, y0, z0, t0, e0, 0.0, 0.0, 0.0);

        int ne = 0, ni = 0;
        aval.GetAvalancheSize(ne, ni);
        if (i % 10 == 0) {
            std::cout << i << "/" << nEvents << "\n"
                    << "    " << ne << " electrons\n";
        }

        #pragma omp critical
        {
            hElectrons.Fill(ne);
            hIons.Fill(ni);
        }

        // // Retrieve the number of collisions of each type.
        // unsigned int nColl[6] = {0};
        // unsigned int nTotal = gas.GetNumberOfElectronCollisions(nColl[0], nColl[1], nColl[2], nColl[3], nColl[4], nColl[5]);
        // nIon += nColl[1];
        // nExc += nColl[4];
    }

    hElectrons.Write();
    hIons.Write();
    // hEnergy.Write();

    const double mean = hElectrons.GetMean();
    const double rms = hElectrons.GetRMS();
    std::cout << " f = " << rms * rms / (mean * mean) << "\n";
    std::printf("    %20u ionisations, %20u excitations\n", nIon, nExc);

    rootOut.Close();

    std::cout << "Done." << std::endl;

    return EXIT_SUCCESS;
}