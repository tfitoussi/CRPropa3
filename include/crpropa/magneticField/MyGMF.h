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
    bool use_disk;
	double B0_D;   // magnetic field scale
	double cos_pitch, sin_pitch, inv_tan_pitch, PHI, cos_PHI;  // pitch angle parameters
	double R_sun; // distance between sun and galactic center
    double R_c; // constant central region

    // halo parameters
    bool use_halo;
	double B0_H;   // magnetic field scale
    double rho_0;
    double rho_1;

    // transition disk - halo
    double z_0; // altitude transition disk halo
    double l_dh; // transition disk halo 

    double transition(const double z) const;

public:
    double logisticalFunction(const double &x, const double &x0, const double &lambda) const;

	MyGMF(
        // disk
        const double B0d, // muG
        const double pitch, // rad
        const double Rc, // kpc
        const double d, // kpc
        // halo
        const double B0h, // muG
        const double rho0, // kpc
        const double rho1, // kpc
        // transition
        const double z0, // kpc
        const double ldh // kpc
    );

	Vector3d getField(const Vector3d& pos) const;

    void setUseDisk(bool use);
	void setUseHalo(bool use);
};

} // namespace crpropa

#endif // CRPROPA_MYGMF_H
