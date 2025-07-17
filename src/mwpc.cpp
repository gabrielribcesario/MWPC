#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <omp.h>

#include <TROOT.h>
#include <TRandom3.h>
#include <TObject.h>
#include <TNtupleD.h>
#include <TFile.h>

#include <Garfield/SolidWire.hh>
#include <Garfield/GeometrySimple.hh>
#include <Garfield/ComponentNeBem3d.hh>
#include <Garfield/Sensor.hh>
#include <Garfield/MediumConductor.hh>
#include <Garfield/MediumMagboltz.hh>
#include <Garfield/DriftLineRKF.hh>
#include <Garfield/TrackHeed.hh>

#include <wires.hpp>

#define KALPHA1_P 16.48 // probability of K_α1
#define KALPHA2_P 8.40 // probability of K_α2
#define KBETA_P 3.38 // probability of K_β
#define KALPHA1_E 5.89881E3 // energy [eV] of K_α1
#define KALPHA2_E 5.88772E3 // energy [eV] of K_α1
#define KBETA_E 6.49051E3 // energy [eV] of K_α1
#define BUFFER_SIZE 2048

using namespace Garfield;

int main(int argc, char **argv) {
    // Random number generator, seed = 42
    TRandom3 rng(42);

    // Number of events to simulate
    const size_t nEvents = 1;

    const char *rootPath = "output.root";
    const char *gasPath = "../gas/1atm/Ar_80_CO2_20.gas";
    const char *ionPath = "IonMobility_Ar+_Ar.txt";

    // Output .root file
    TFile rootOut(rootPath, "create", "", 505);

    // Load electron and ion transport parameters
    MediumMagboltz gas;
    if (!gas.LoadGasFile(gasPath) || !gas.LoadIonMobility(ionPath)) { return EXIT_FAILURE; }
    gas.EnablePenningTransfer();

    // Collimation radius [cm]
    const double collimation = 0.05;
    // Cathode wire spacing [cm]
    const double catG = 0.2;
    // Cathode wire radius [cm]
    const double catR = 0.005;
    // Cathode wire potential [V]
    const double catV = 0.;
    // Anode wire spacing [cm]
    const double anoG = 0.4;
    // Anode wire diameter [cm]
    const double anoR = 0.002;
    // Anode wire potential [V]
    const double anoV = 2100.;
    // Cathodic-Anodic planes distance [cm], anodic plane at z = 0
    const double catZ = 0.5;
    // Wire length [cm], at least XX gaps in both the x and y directions
    const double wLen = std::ceil(anoG > catG ? 8.0 / anoG : 8.0 / catG);

    // Cathode wires plane, upper
    std::vector<SolidWire> catPlaneU;
    if(!SetWires(catPlaneU, static_cast<int>(std::floor(wLen / catG * 0.5)),
                 0.0, catG, catZ, catR, wLen, 1.0, 0.0)) { return EXIT_FAILURE; }
    // Cathode wires plane, lower
    std::vector<SolidWire> catPlaneL;
    if (!SetWires(catPlaneL, static_cast<int>(std::floor(wLen / catG * 0.5)),
                  0.0, catG, -catZ, catR, wLen, 1.0, 0.0)) { return EXIT_FAILURE; }
    // // Anode wires plane
    std::vector<SolidWire> anoPlane;
    if (!SetWires(anoPlane, static_cast<int>(std::floor(wLen / anoG * 0.5)), 
                  anoG, 0.0, 0.0, anoR, wLen, 0.0, 1.0)) { return EXIT_FAILURE; }

    MediumConductor metal;

    // Set the geometry, electrodes and electric potential
    GeometrySimple geo;
    geo.SetMedium(&gas);
    AddWires(geo, metal, catPlaneU, catV, "catU", nullptr);
    AddWires(geo, metal, catPlaneL, catV, "catL", nullptr);
    AddWires(geo, metal, anoPlane, anoV, "ano", nullptr);

    // Electric field calculation
    ComponentNeBem3d cmp;
    cmp.SetGeometry(&geo);
    cmp.SetTargetElementSize(0.001);
    cmp.SetMinMaxNumberOfElements(1, 5); // fine-tune element size and number later
    cmp.UseLUInversion();
    cmp.SetNumberOfThreads(8);
    if (!cmp.Initialise()) { return EXIT_FAILURE; }

    // Simulation boundaries
    const double xlim = wLen * 0.5;
    const double ylim = wLen * 0.5;
    const double zlim = catZ * 1.25;
    // Interface between the transport classes and the component
    Sensor sensor;
    sensor.AddComponent(&cmp);
    sensor.SetArea(-xlim, -ylim, -zlim,
                    xlim,  ylim,  zlim);

    // DriftLineRKF maximum step size
    const double stepRKF = 0.01;
    std::printf("DriftLineRKF max step size [cm]: %#g\n", stepRKF);
    // DriftLineRKF integration accuracy
    const double epsRKF = 1.0E-8;
    std::printf("DriftLineRKF integration accuracy: %#g\n", epsRKF);
    // Charge transport
    DriftLineRKF driftRKF(&sensor);
    driftRKF.SetMaximumStepSize(stepRKF);
    driftRKF.SetIntegrationAccuracy(epsRKF);

    // Primary ionizations calculations
    TrackHeed track(&sensor);

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
    double dx = 0.0, dy = 0.0, dz = 0.0;
    // Photon track tabular data
    TNtupleD trackTree("Photons", "Photon track data", "x0:y0:z0:dx:dy:dz:e0:np");

    for (unsigned int event = 1; event <= nEvents; ++event) {  
        int ntries = 0;
        // Discards the event if the photon does not interact with the detector
        while (!np && ++ntries < static_cast<int>(KBETA_E)) {
            rng.Circle(x0, y0, collimation); // replicate collimation
            rng.Sphere(dx, dy, dz, 1.0); // normalized direction
            dx = std::signbit(dx * x0) ? -dx : dx; // same sign (entry point on plane A, vector should lie on plane A)
            dy = std::signbit(dy * y0) ? -dy : dy; // same sign (entry point on plane A, vector should lie on plane A)
            dz = std::signbit(dz) ? -dz : dz; // always pointing down
            z0 = 1.1 * (catZ + catR * 0.5); // fixed z0
            e0 = rng.Uniform(0.0, KALPHA1_P + KALPHA2_P + KBETA_P); // sample x-ray energy
            e0 = e0 < KBETA_P ? KBETA_E : e0 < KBETA_P + KALPHA1_P ? KALPHA1_E : KALPHA2_E; 
            track.TransportPhoton(x0, y0, z0, 0.0, e0, dx, dy, dz, np);
        }
        if (ntries >= static_cast<int>(KBETA_E)) { 
            std::cerr << "Could not find an initial position with an ionisable medium." << std::endl;
            return EXIT_FAILURE; 
        }

        trackTree.Fill(x0, y0, z0, dx, dy, dz, e0, static_cast<double>(np));
        trackTree.Write(nullptr, TObject::kWriteDelete);

        std::cout << "Event #" << event << std::endl;
        std::cout << "|   Photon energy [eV]: " << e0 << std::endl;
        std::cout << "|   Photon conversion: " << np << "electron-ion pairs" << std::endl;
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