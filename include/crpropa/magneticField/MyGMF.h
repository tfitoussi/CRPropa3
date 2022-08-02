#ifndef CRPROPA_MYGMF_H
#define CRPROPA_MYGMF_H

#include "crpropa/magneticField/MagneticField.h"

namespace crpropa {

/**
 @class MyGMF
 @brief Implements my own galactic magnetic field model
*/

class MyGMF: public MagneticField {
private:
	// disk parameters
	double B_0d;   // magnetic field scale
    bool use_disk;
	double cos_pitch, sin_pitch, tan_pitch;  // pitch angle parameters
	double R_sun; // distance between sun and galactic center
    double Phi_0; // disk phase

    // halo parameters
    bool use_halo;
	double B_0h;   // magnetic field scale
    double rho_0;
    double rho_1;

    // transition disk - halo
    double R_gc; // radius central region
    double l_gc; // transition central region / disk
    double R_b; // scale exponential decay disk at Earth
    double z_0; // altitude transition disk halo
    double l_dh; // transition disk halo 

    double transition(const double r, const double z) const;

public:
    double logisticalFunction(const double &x, const double &x0, const double &lambda) const;

	MyGMF(
        // disk
        const double B0d, // muG
        const double pitch, // rad
        const double phi0, // rad
        // halo
        const double B0h, // muG
        const double rho0, // kpc
        const double rho1, // kpc
        // transition
        const double Rgc, // kpc
        const double lgc, // kpc
        const double Rb, // kpc
        const double z0, // kpc
        const double ldh // kpc
    );

	Vector3d getField(const Vector3d& pos) const;

    void setUseDisk(bool use);
	void setUseHalo(bool use);
};

} // namespace crpropa

#endif // CRPROPA_MYGMF_H
