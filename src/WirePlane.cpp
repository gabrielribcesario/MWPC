#include "WirePlane.hpp"

#include <iostream>
#include <string>

#include <Garfield/GeometrySimple.hh>
#include <Garfield/SolidWire.hh>
#include <Garfield/MediumMagboltz.hh>
#include <Garfield/MediumConductor.hh>

bool WirePlane::SetWirePlane(const unsigned long ngaps,
                             const double z0, const double gap, const double radius, const double potential, 
                             const std::string &orientation,
                             Garfield::GeometrySimple *geo, Garfield::MediumConductor *metal) {
    if (!geo) { 
        std::cerr << "WirePlane::SetWirePlane: Garfield::GeometrySimple* is nullptr\n";
        return false;
    }
    if (!metal) {
        std::cerr << "WirePlane::SetWirePlane: Garfield::MediumConducto* is nullptr\n";
        return false;
    }

    double dx, dy, dz;
    if (orientation == "x") { dx = 1.0; dy = 0.0; dz = 0.0; }
    else if (orientation == "y") { dx = 0.0; dy = 1.0; dz = 0.0; }
    else {
        std::cerr << "WirePlane::SetWirePlane Error: Unknown wire orientation, options are 'x' and 'y'\n";
        return false;
    }

    m_orientation = orientation.c_str();
    m_ngaps = ngaps;
    m_z0 = z0;
    m_gap = gap;
    m_radius = radius;
    m_potential = potential;

    m_wires.reserve(m_ngaps + 1);

    const double length = m_ngaps * m_gap;

    // Add the wires to the list of primitives
    for (unsigned long i = 0; i <= m_ngaps; ++i) {
        double x0, y0;
        if (dx == 1.0) {
            x0 = 0.0;
            y0 = (m_ngaps * -0.5 + i) * m_gap;
        }
        else {
            x0 = (m_ngaps * -0.5 + i) * m_gap;
            y0 = 0.0;
        }

        m_wires.emplace_back(x0, y0, m_z0, m_radius, length * 0.5, dx, dy, dz);
        m_wires.back().SetBoundaryPotential(m_potential);
        m_wires.back().SetColour(m_colour);
        geo->AddSolid(&m_wires.back(), metal);
    }

    return true;
}

void WirePlane::SetColour(const int colour) {
    m_colour = colour;
    for (auto &wire : m_wires) {
        wire.SetColour(m_colour);
    }
}