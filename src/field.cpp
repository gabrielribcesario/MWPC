#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>

#include <TROOT.h>
#include <TApplication.h>

#include <Garfield/SolidWire.hh>
#include <Garfield/GeometrySimple.hh>
#include <Garfield/ComponentNeBem3d.hh>
#include <Garfield/MediumConductor.hh>
#include <Garfield/MediumMagboltz.hh>
#include <Garfield/ViewField.hh>
#include <Garfield/ViewGeometry.hh>

#include <wires.hpp>

#define BUFFER_SIZE 2048

using namespace Garfield;

int main(int argc, char **argv) {
    // Cathode wire spacing [cm].
    const double catG = 0.2;
    // Cathode wire radius [cm].
    const double catR = 0.005;
    // Cathode wire potential [V].
    const double catV = 0.;
    // Anode wire spacing [cm].
    const double anoG = 0.4;
    // Anode wire diameter [cm].
    const double anoR = 0.002;
    // Anode wire potential [V].
    const double anoV = 2100.;
    // Cathodic-Anodic planes distance [cm], anodic plane at z = 0
    const double catZ = 0.5;
    // Wire length [cm].
    const double wLen = 3.0;

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
    MediumMagboltz gas;

    int colour = kBlack;
    // Set the geometry, electrodes and electric potential
    GeometrySimple geo;
    geo.SetMedium(&gas);
    AddWires(geo, metal, catPlaneU, catV, nullptr, &colour);
    AddWires(geo, metal, catPlaneL, catV, nullptr, &colour);
    colour = kYellow + 1;
    AddWires(geo, metal, anoPlane, anoV, nullptr, &colour);

    ComponentNeBem3d cmp;
    cmp.SetGeometry(&geo);
    cmp.SetTargetElementSize(0.001);
    cmp.SetMinMaxNumberOfElements(1, 5); // fine-tune element size and number later
    cmp.UseLUInversion();
    cmp.SetNumberOfThreads(8);
    cmp.Initialise();

    TApplication app("app", &argc, argv);

    const double xlim = catG * 2.;
    const double ylim = anoG * 2.;
    const double zlim = catZ * 1.1;

    // XZ view of the electric potential field
    ViewField fieldViewXZ(&cmp);
    fieldViewXZ.SetArea(-xlim, -zlim,
                         xlim,  zlim);
    fieldViewXZ.SetPlaneXZ();
    fieldViewXZ.PlotContour();

    // YZ view of the electric potential field
    ViewField fieldViewYZ(&cmp);
    fieldViewYZ.SetArea(-ylim, -zlim,
                         ylim,  zlim);
    fieldViewYZ.SetPlaneYZ();
    fieldViewYZ.PlotContour();

    ViewGeometry geoView(&geo);
    geoView.SetArea(-xlim, -ylim, -zlim, 
                     xlim,  ylim,  zlim);
    geoView.Plot3d();

    std::cout << "Done." << std::endl;

    app.Run(true);

    return EXIT_SUCCESS;
}