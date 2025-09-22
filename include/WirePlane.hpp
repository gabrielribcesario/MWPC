#ifndef WIREPLANE_HPP
#define WIREPLANE_HPP

#include <iostream>
#include <string>
#include <memory>
#include <vector>

namespace Garfield {
    class GeometrySimple;
    class SolidWire;
    class MediumConductor;
}

// Encapuslation of the geomery of a plane of parallel solid wires
class WirePlane {
    public:
    WirePlane() = default;

    WirePlane(const unsigned long ngaps, 
              const double z0, const double gap, const double radius, const double potential, 
              const std::string &orientation, 
              Garfield::GeometrySimple *geo, Garfield::MediumConductor *metal) {
        SetWirePlane(ngaps, z0, gap, radius, potential, orientation, geo, metal);
    }

    ~WirePlane() = default;

    bool SetWirePlane(const unsigned long ngaps, 
                      const double z0, const double gap, const double radius, const double potential, 
                      const std::string &orientation,
                      Garfield::GeometrySimple *geo, Garfield::MediumConductor *metal);

    // Set the plotting colour of the wires
    void SetColour(const int colour);

    private:
    const char *m_className = "WirePlane";

    const char *m_orientation;

    std::vector<Garfield::SolidWire> m_wires;

    // Number of wire gaps
    unsigned long m_ngaps;

    // Centroid of the wire plane along the z-axis [cm]
    double m_z0;

    // Wire spacing [cm]
    double m_gap;

    // Wire radius [cm]
    double m_radius;

    // Wire electric potential [V]
    double m_potential;

    // Wire plotting colour
    int m_colour = -1;
};

#endif