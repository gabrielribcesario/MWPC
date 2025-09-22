#include <cstdlib>
#include <vector>
#include <iostream>

#include <TApplication.h>
#include <TColor.h>

#include <Garfield/MediumConductor.hh>
#include <Garfield/MediumMagboltz.hh>
#include <Garfield/GeometrySimple.hh>
#include <Garfield/SolidBox.hh>
#include <Garfield/SolidWire.hh>
#include <Garfield/ComponentNeBem3d.hh>
#include <Garfield/MediumConductor.hh>
#include <Garfield/MediumMagboltz.hh>
#include <Garfield/ViewField.hh>
#include <Garfield/ViewGeometry.hh>

#include "WirePlane.hpp"

using namespace Garfield;

int main(int argc, char **argv) {
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
    // Cathode plotting colour
    constexpr int cathodeColour = kBlack;

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
    // Anode plotting colour
    constexpr int anodeColour = kYellow + 1;

    // Cathode wire length [cm]
    constexpr double cathodeLength = anodeNGaps * anodeGap;
    // Anode wire length [cm]
    constexpr double anodeLength = cathodeNGaps * cathodeGap;

    // Wire material
    MediumConductor metal;

    // Medium
    MediumMagboltz gas;

    // Detector geometry
    GeometrySimple geo;
    geo.SetMedium(&gas);

    WirePlane cathodeUpperPlane(cathodeNGaps, planeGap, cathodeGap, 
                                cathodeRadius, cathodePotential, cathodeOrientation, 
                                &geo, &metal);
    cathodeUpperPlane.SetColour(cathodeColour);

    WirePlane anodePlane(anodeNGaps, 0.0, anodeGap, 
                         anodeRadius, anodePotential, anodeOrientation, 
                         &geo, &metal);
    anodePlane.SetColour(anodeColour);

    WirePlane cathodeLowerPlane(cathodeNGaps, -planeGap, cathodeGap, 
                                cathodeRadius, cathodePotential, cathodeOrientation, 
                                &geo, &metal);
    cathodeLowerPlane.SetColour(cathodeColour);

    ComponentNeBem3d cmp;
    cmp.SetGeometry(&geo);
    cmp.SetTargetElementSize(0.1);
    cmp.UseLUInversion();
    cmp.SetNumberOfThreads(36);
    // cmp.SetReuseModel(); // reuse mesh and model, new boundary conditions
    cmp.Initialise();

    TApplication app("app", &argc, argv);

    constexpr double xlim = cathodeGap * cathodeNGaps * 0.5;
    constexpr double ylim = anodeGap * anodeNGaps * 0.5;
    constexpr double zlim = planeGap * 1.1;

    // XZ view of the electric potential field
    ViewField fieldViewXZ(&cmp);
    fieldViewXZ.SetArea(-xlim, -zlim,
                         xlim,  zlim);
    fieldViewXZ.SetPlaneXZ();
    fieldViewXZ.PlotContour("v");
    fieldViewXZ.GetCanvas()->SetTitle("Electric potential contour, X-Z view");
    fieldViewXZ.GetCanvas()->Modified();
    fieldViewXZ.GetCanvas()->Update();

    // YZ view of the electric potential field
    ViewField fieldViewYZ(&cmp);
    fieldViewYZ.SetArea(-ylim, -zlim,
                         ylim,  zlim);
    fieldViewYZ.SetPlaneYZ();
    fieldViewYZ.PlotContour("v");
    fieldViewYZ.GetCanvas()->SetTitle("Electric potential contour, Y-Z view");
    fieldViewYZ.GetCanvas()->Modified();
    fieldViewYZ.GetCanvas()->Update();

    // 3D view of the geometry
    ViewGeometry geoView(&geo);
    geoView.SetArea(-xlim, -ylim, -zlim, 
                     xlim,  ylim,  zlim);
    geoView.Plot3d();
    geoView.GetCanvas()->SetTitle("3D view of the geometry");
    geoView.GetCanvas()->Modified();
    geoView.GetCanvas()->Update();

    std::cout << "Done." << std::endl;

    app.Run(true);

    return EXIT_SUCCESS;
}