#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TNtupleD.h>

#include <Garfield/MediumConductor.hh>
#include <Garfield/MediumMagboltz.hh>
#include <Garfield/GeometrySimple.hh>
#include <Garfield/SolidWire.hh>
#include <Garfield/ComponentNeBem3d.hh>
#include <Garfield/Sensor.hh>
#include <Garfield/DriftLineRKF.hh>
#include <Garfield/TrackHeed.hh>
#include <Garfield/Random.hh>

#include "WirePlane.hpp"

static constexpr double fe55_kalpha_prob = 24.88; // relative probability of k-alpha x-ray
static constexpr double fe55_kbeta_prob = 3.38; // relative probability of k-beta x-ray

static constexpr double fe55_kalpha_eV = 5.9e03; // [eV]
static constexpr double fe55_kbeta_eV = 6.49e03; // [eV]

using namespace Garfield;

int main(int argc, char **argv) {
    const char *rootPath = "output.root";
    const char *gasPath = "../data/gas/1atm/Ar_80_CO2_20.gas";
    const char *ionPath = "IonMobility_Ar+_Ar.txt";

    // Number of events to simulate
    const size_t nEvents = 10;

    // DriftLineRKF maximum step size
    const double stepRKF = 0.01;
    // DriftLineRKF integration accuracy
    const double epsRKF = 1.0E-7;

    // Cathode-Anode planes spacing [cm]
    constexpr double planeGap = 0.5;

    // Wire orientation
    const char *cathodeOrientation = "x";
    // Number of cathode wire gaps
    constexpr size_t cathodeNGaps = 8;
    // Cathode wire spacing [cm]
    constexpr double cathodeGap = 0.2;
    // Cathode wire radius [cm]
    constexpr double cathodeRadius = 0.0025;
    // Cathode wire potential [V]
    constexpr double cathodePotential = 0.;

    // Wire orientation
    const char *anodeOrientation = "y";
    // Number of anode wire gaps
    constexpr size_t anodeNGaps = 4;
    // Anode wire spacing [cm]
    constexpr double anodeGap = 0.4;
    // Anode wire radius [cm]
    constexpr double anodeRadius = 0.001;
    // Anode wire potential [V]
    constexpr double anodePotential = 2100.;

    // Cathode wire length [cm]
    constexpr double catWireLen = anodeNGaps * anodeGap;
    // Anode wire length [cm]
    constexpr double anoWireLen = cathodeNGaps * cathodeGap;

    // Wire material
    MediumConductor metal;

    // Medium
    MediumMagboltz gas;
    if (!gas.LoadGasFile(gasPath) || !gas.LoadIonMobility(ionPath)) { return EXIT_FAILURE; }
    gas.EnablePenningTransfer();

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
    if (!cmp.Initialise()) { return EXIT_FAILURE; }

    // std::cout   << "Geometry" << std::endl
    //             << "|   Anode wires" << std::endl
    //             << "|   |   Count: " << anoPlane.size() << " wires" << std::endl
    //             << "|   |   Spacing: " << anoG << "[cm]" << std::endl
    //             << "|   |   Diameter: " << anoR * 2.0 << "[cm]" << std::endl
    //             << "|   |   Electric potential: " << anoV << "[V]" << std::endl
    //             << "|   Cathode wires" << std::endl
    //             << "|   |   Count: " << catPlaneU.size() * 2 << " wires" << std::endl
    //             << "|   |   Spacing: " << catG << "[cm]" << std::endl
    //             << "|   |   Diameter: " << catR * 2.0 << "[cm]" << std::endl
    //             << "|   |   Electric potential: " << catV << "[V]" << std::endl;

    // Bounding box
    constexpr double xlim = cathodeGap * cathodeNGaps * 0.5;
    constexpr double ylim = anodeGap * anodeNGaps * 0.5;
    constexpr double zlim = planeGap * 1.1;

    // Interface between the transport classes and the component
    Sensor sensor;
    sensor.AddComponent(&cmp);
    sensor.SetArea(-xlim, -ylim, -zlim,
                    xlim,  ylim,  zlim);

    // Charge transport
    DriftLineRKF driftRKF(&sensor);
    driftRKF.SetMaximumStepSize(stepRKF);
    driftRKF.SetIntegrationAccuracy(epsRKF);

    std::cout   << "DriftLineRKF\n"
                << "|   Integration accuracy: " << epsRKF << "\n"
                << "|   Maximum step size: " << stepRKF << std::endl;

    // Primary ionizations calculations
    TrackHeed track(&sensor);

    // Output .root file
    TFile rootOut(rootPath, "recreate", "", 505);

    // Drift line status flag
    int status = 0;
    // Drift line gain
    double gain = 0.0;
    // Drift line loss
    double loss = 0.0;
    // # of electrons at the end of the drift line
    double ne = 0.0;
    // # of ions at the end of the ion tail
    double ni = 0.0;
    // Drift line tabular data
    TNtupleD driftTree("DriftLines", "Drift line data", "status:gain:loss:ne:ni");

    // Initial point [cm]
    double x0 = 0.0, y0 = 0.0, z0 = 0.0;
    // Photon energy [eV]
    double e0 = 0.0;
    // # of electrons (photon conversion)
    int np = 0;
    // Photon direction
    double dx = 0.0, dy = 0.0, dz = -1.0;
    // Photon track tabular data
    TNtupleD trackTree("Photons", "Photon track data", "x0:y0:z0:dx:dy:dz:e0:np");

    for (unsigned int event = 1; event <= nEvents; ++event) {  
        np = 0;
        int ntries = 0;
        // Discards the event if the photon does not interact with the detector
        while (!np && ++ntries < static_cast<int>(fe55_kbeta_eV)) {
            double rms = 0.05;

            // Replicate the collimation process
            x0 = RndmGaussian(0.0, rms);;
            y0 = RndmGaussian(0.0, rms);;
            z0 = planeGap - cathodeRadius - 1e-3; // fixed z0

            double angle = 45.;

            // do { // Uniform spherical distribution with truncated incidence angle
            //     rng.Sphere(dx, dy, dz, 1.0); 
            // } while (std::atan2(std::sqrt( dx * dx + dy * dy ), dz) > angle * 3.141592653589793 / 180.0);

            dx = dx * x0 > 0.0 ? -dx : dx; // same sign (entry point on plane A, vector should lie on plane A)
            dy = dy * y0 > 0.0 ? -dy : dy; // same sign (entry point on plane A, vector should lie on plane A)
            dz = dz > 0.0 ? dz : -dz; // always pointing down

            // Sample x-ray energy
            e0 = RndmUniform() < (fe55_kalpha_prob / (fe55_kalpha_prob + fe55_kbeta_prob)) ? fe55_kalpha_eV : fe55_kbeta_eV;
            track.TransportPhoton(x0, y0, z0, 0.0, e0, dx, dy, dz, np);
        }
        if (ntries >= static_cast<int>(fe55_kbeta_eV)) { 
            std::cerr << "Could not find an initial position with an ionisable medium\n";
            return EXIT_FAILURE; 
        }

        trackTree.Fill(x0, y0, z0, dx, dy, dz, e0, static_cast<double>(np));
        trackTree.Write(nullptr, TObject::kWriteDelete);

        std::cout << "Event #" << event << std::endl;
        std::cout << "|   Photon energy [eV]: " << e0 << std::endl;
        std::cout << "|   Photon conversion: " << np << " electron-ion pairs" << std::endl;
        std::cout << "|   (x0, y0, z0): (" << x0 << ", " << y0 << ", " << z0 << ")" << std::endl;
        std::cout << "|   (dx, dy, dz): (" << dx << ", " << dy << ", " << dz << ")" << std::endl;

        // Loop over the primary electrons
        for (const auto &cluster : track.GetClusters()) {
            #pragma omp parallel for firstprivate(driftRKF) private(status,gain,loss,ne,ni)
            for (const auto &electron : cluster.electrons) {
                driftRKF.DriftElectron(electron.x, electron.y, electron.z, electron.t);

                // Get drift line endpoint status
                double x1 = 0.0, y1 = 0.0, z1 = 0.0, t1 = 0.0;
                driftRKF.GetEndPoint(x1, y1, z1, t1, status);

                // Integrate the drift line outside the critical region
                gain = driftRKF.GetGain();
                loss = driftRKF.GetLoss();

                // Get avalanche size
                driftRKF.GetAvalancheSize(ne, ni);

                #pragma omp critical 
                {
                    driftTree.Fill(static_cast<double>(status), gain, loss, ne, ni);
                }
            }
        }
        driftTree.Write(nullptr, TObject::kWriteDelete);
    }
    rootOut.Close();

    std::cout << "Done." << std::endl;

    return EXIT_SUCCESS;
}