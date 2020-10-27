#ifndef CRPROPA_FMO20FIELD_H
#define CRPROPA_FMO20FIELD_H

#include "crpropa/magneticField/MagneticField.h"

namespace crpropa {

class FMO20Field: public MagneticField {
private:
    // Disk parameters
    bool use_disk;
    double B0_disk;
    double Rc_disk;
    double d_disk;
    double pitch_disk, sin_pitch_disk, cos_pitch_disk;
    double PHI, cos_PHI;

    // Halo parameters
    bool use_halo;
    double B0_halo;
    double rho0_halo;
    double rho1_halo;

    // transition disk - halo
    double z0;
    double k_LF;

    // Global parameter
    double R_sun;
    

public:
    FMO20Field();

	Vector3d getField(const Vector3d& pos) const;
	Vector3d getDiskField(const Vector3d& pos) const;
	Vector3d getHaloField(const Vector3d& pos) const;

    double logisticalFunction(const double& z, const double& z0, const double& k) const;

    void setB0Disk(double B);
    void setRcDisk(double Rc);
    void setDDisk(double d);
    void setPitchAngle(double pitch_angle);
    void setParams();

    void setB0Halo(double B);
    void setRho0Halo(double rho);
    void setRho1Halo(double rho);

    void setZ0(double z);
    void setKLF(double k);

	void setUseDisk(bool use);
	void setUseHalo(bool use);
	bool isUsingDisk();
	bool isUsingHalo();
};

} // namespace crpropa

#endif // CRPROPA_JF12FIELDSOLENOIDAL_H
