#ifndef NEBEMWIRES_HPP
#define NEBEMWIRES_HPP

#include <cstdlib>
#include <iostream>
#include <vector>

#include <Garfield/SolidWire.hh>
#include <Garfield/GeometrySimple.hh>
#include <Garfield/MediumConductor.hh>

using namespace Garfield;

// Add the wires to the geometry
bool SetWires(std::vector<SolidWire> &wires, const int ngap,
              const double x0, const double y0, const double z0, 
              const double radius, const double length,
              const double dx, const double dy) {
    if ((x0 != 0.0) == (y0 != 0.0)) { // XOR: either x0 or y0, neither equal to zero
        std::cerr   << "SetWires() failure: Either x0 or y0 should be set." << std::endl 
                    << "(x0, y0) = (" << x0 << ", " << y0 << ")" << std::endl;
        return false; 
    }
    wires.reserve(ngap * 2 + 1);
    wires.emplace_back(0.0, 0.0, z0,
                       radius, length * 0.5, 
                       dx, dy, 0.0); // Wire #0 lies on (0,0)
    for (int i = 1; i <= ngap; ++i) { 
        wires.emplace_back(x0 * i, y0 * i, z0,
                           radius, length * 0.5, 
                           dx, dy, 0.0); // Odd wires lie on (x,y) > (0,0)
        wires.emplace_back(-x0 * i, -y0 * i, z0,
                           radius, length * 0.5, 
                           dx, dy, 0.0); // Even wires lie on (x,y) < (0,0)
    }
    return true;
}

// Add wires, set electrode label and electric potential
bool AddWires(GeometrySimple &geo, MediumConductor &mat, std::vector<SolidWire> &wires, 
              const double potV, const char *label, const int *colour) {
    int w = -1;
    for (auto &wire : wires) {
        if (label) {
            // std::string buffer = (std::string) label + std::to_string(++w);
            std::string buffer = (std::string) "wire" + std::to_string(wire.GetId());
            wire.SetLabel(buffer); 
        }
        if (colour) { wire.SetColour(*colour); }
        wire.SetBoundaryPotential(potV);
        geo.AddSolid(&wire, &mat);
    }
    return true;
}

#endif