#include <cstdlib>
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
    // Anode-Cathode wires gap [cm];
    const double wGap = 0.5;
    // Wire length [cm].
    const double wLen = 10.;

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

    // Cathode wires plane, upper
    SolidWire catWireU(0., 0., wGap, 
                       catR, wLen * 0.5, 
                       1., 0., 0.);
    catWireU.SetBoundaryPotential(catV);

    // Cathode wires plane, lower
    SolidWire catWireL(0., 0., -wGap, 
                       catR, wLen * 0.5, 
                       1., 0., 0.);
    catWireL.SetBoundaryPotential(catV);

    // Anode wires plane
    SolidWire anoWire(0., 0., 0., 
                      anoR, wLen * 0.5, 
                      0., 1., 0.);
    anoWire.SetBoundaryPotential(anoV);

    MediumConductor metal;
    MediumMagboltz gas;

    GeometrySimple geo;
    geo.SetMedium(&gas);
    geo.AddSolid(&catWireU, &metal);
    geo.AddSolid(&catWireL, &metal);
    geo.AddSolid(&anoWire, &metal);

    ComponentNeBem3d cmp;
    cmp.SetGeometry(&geo);
    cmp.SetPeriodicityX(anoG); // no signal calculation --> use periodicity
    cmp.SetPeriodicityY(catG); // no signal calculation --> use periodicity
    cmp.SetTargetElementSize(0.001);
    cmp.SetMinMaxNumberOfElements(3, 15); // fine-tune element size and number later
    cmp.UseLUInversion();
    cmp.SetNumberOfThreads(8);
    cmp.Initialise();

    TApplication app("app", &argc, argv);
    
    // XZ view of the electric potential field
    ViewField fieldViewXZ(&cmp);
    fieldViewXZ.SetArea(-catG * 2., -wGap,
                         catG * 2.,  wGap);
    fieldViewXZ.SetPlaneXZ();
    fieldViewXZ.PlotContour();

    // YZ view of the electric potential field
    ViewField fieldViewYZ(&cmp);
    fieldViewYZ.SetArea(-catG * 2., -wGap,
                         catG * 2.,  wGap);
    fieldViewYZ.SetPlaneYZ();
    fieldViewYZ.PlotContour();

    // ViewGeometry geoView(&geo);
    // geoView.SetArea(-catG * 2., -anoG * 2., -wGap * 1.2, 
    //                  catG * 2.,  anoG * 2.,  wGap * 1.2);
    // geoView.Plot3d();

    std::cout << "Done." << std::endl;

    app.Run(true);

    return EXIT_SUCCESS;
}