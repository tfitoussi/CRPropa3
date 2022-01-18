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
	double B_0;   // magnetic field scale

	// disk parameters
    bool use_disk;
	double cos_pitch, sin_pitch, PHI, cos_PHI;  // pitch angle parameters
	double d;     // distance to first field reversal
	double R_sun; // distance between sun and galactic center

    // halo parameters
    bool use_halo;
    double rho_0;
    double rho_1;

    // transition disk - halo
	double R_b;   // radius of central region
	double z_0;   // vertical thickness in the galactic disk
    double lambda_gc, lambda_r, lambda_z;

    double transition(const double r, const double z) const;

public:
    double logisticalFunction(const double &x, const double &x0, const double &lambda) const;

	MyGMF(const double B, const double R, const double z, const double D, const double pitch,
            const double rho0, const double rho1, const double lambda_r, const double lambda_z);

	Vector3d getField(const Vector3d& pos) const;

    void setUseDisk(bool use);
	void setUseHalo(bool use);
};

} // namespace crpropa

#endif // CRPROPA_MYGMF_H
