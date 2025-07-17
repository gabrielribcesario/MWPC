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
    const double anoV = -2100.;

    // Cathodic-Anodic planes distance [cm];
    const double wGap = 0.5;
    // Wire length [cm].
    const double wLen = 3.0;
    // const double wLen = std::ceil(anoG > catG ? 8.0 / anoG : 8.0 / catG); // at least 8 gaps in the x and y directions

    // Cathode wires plane, upper
    std::vector<SolidWire> catPlaneU;
    catPlaneU.emplace_back(0., 0., wGap,
                           catR, wLen * 0.5, 
                           1., 0., 0.);

    // Cathode wires plane, lower
    std::vector<SolidWire> catPlaneL;
    catPlaneU.emplace_back(0., 0., -wGap,
                           catR, wLen * 0.5, 
                           1., 0., 0.);

    // // Anode wires plane
    std::vector<SolidWire> anoPlane;
    anoPlane.emplace_back(0., 0., 0.,
                          anoR, wLen * 0.5, 
                          0., 1., 0.);

    for (int i = 1; i <= (int) std::floor(wLen / catG * 0.5); ++i) {
        catPlaneU.emplace_back(0., catG * i, wGap,
                               catR, wLen * 0.5, 
                               1., 0., 0.);
        catPlaneU.emplace_back(0., -catG * i, wGap,
                               catR, wLen * 0.5, 
                               1., 0., 0.);
        catPlaneL.emplace_back(0., catG * i, -wGap,
                               catR, wLen * 0.5, 
                               1., 0., 0.);
        catPlaneL.emplace_back(0., -catG * i, -wGap,
                               catR, wLen * 0.5, 
                               1., 0., 0.);
    }
    for (int i = 1; i <= (int) std::floor(wLen / anoG * 0.5); ++i) {
        anoPlane.emplace_back(anoG * i, 0., 0.,
                              anoR, wLen * 0.5, 
                              0., 1., 0.);
        anoPlane.emplace_back(-anoG * i, 0., 0.,
                              anoR, wLen * 0.5, 
                              0., 1., 0.);
    }

    MediumConductor metal;
    MediumMagboltz gas;

    GeometrySimple geo;
    geo.SetMedium(&gas);
    for (auto &wire : catPlaneU) {
        wire.SetBoundaryPotential(catV);
        wire.SetColour(kBlack);
        geo.AddSolid(&wire, &metal);
    }
    for (auto &wire : catPlaneL) {
        wire.SetBoundaryPotential(catV);
        wire.SetColour(kBlack);
        geo.AddSolid(&wire, &metal);
    }
    for (auto &wire : anoPlane) {
        wire.SetBoundaryPotential(anoV);
        wire.SetColour(kYellow + 1);
        geo.AddSolid(&wire, &metal);
    }

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
    const double zlim = wGap * 1.1;

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