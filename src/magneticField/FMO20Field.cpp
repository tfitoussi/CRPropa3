#include "crpropa/magneticField/FMO20Field.h"
#include "crpropa/Units.h"
#include <algorithm>
#include <stdio.h>

namespace crpropa {
using namespace std;

FMO20Field::FMO20Field() { 
    // Global parameters
	R_sun = 8.5 * kpc;

    // Disk parameters
    setUseDisk(true); 
    setB0Disk(1*muG);
    setRcDisk(15*kpc);
    setDDisk(-1*kpc);
    setPitchAngle(-12*deg);

    // Halo parameters
    setUseHalo(true);
    setB0Halo(0.95*muG);
    //setRho0Halo(1*kpc);
    setRho1Halo(1*kpc);

    // transition disk - halo
    setZ0(1*kpc);
    setKLF(50./kpc);
}

Vector3d FMO20Field::getField(const Vector3d& pos) const {
    Vector3d B(0.);
    double LF = logisticalFunction(pos.z, z0, k_LF);
    if (use_disk) B += getDiskField(pos) * (1.-LF);
    if (use_halo) B += getHaloField(pos) * LF;
    return B;
}

Vector3d FMO20Field::getDiskField(const Vector3d& pos) const {
	double r = sqrt(pos.x * pos.x + pos.y * pos.y);  // in-plane radius
    // Define using Pshirkov et al. 2011 bisymmetrical disk field (BSS)
    Vector3d B(0.);

    // PT11 paper has B_theta = B * cos(p) but this seems because they define azimuth clockwise, while we have anticlockwise.
    // see Tinyakov 2002 APh 18,165: "local field points to l=90+p" so p=-5 deg gives l=85 and hence clockwise from above.
    // so to get local B clockwise in our system, need minus (like Sun et al).
    // Ps base their system on Han and Qiao 1994 A&A 288,759 which has a diagram with azimuth clockwise, hence confirmed.

    // PT11 paper define Earth position at (+8.5, 0, 0) kpc; but usual convention is (-8.5, 0, 0)
    // thus we have to rotate our position by 180 degree in azimuth
    double theta = M_PI - pos.getPhi();  // azimuth angle theta: PT11 paper uses opposite convention for azimuth
    // the following is equivalent to sin(pi - phi) and cos(pi - phi) which is computationally slower
    double cos_theta = - pos.x / r;
    double sin_theta = pos.y / r;

    // After some geometry calculations (on whiteboard) one finds:
    // Bx = +cos(theta) * B_r - sin(theta) * B_{theta}
    // By = -sin(theta) * B_r - cos(theta) * B_{theta}
    // Use from paper: B_theta = B * cos(pitch)	and B_r = B * sin(pitch)
    B.x = sin_pitch_disk * cos_theta - cos_pitch_disk * sin_theta;
    B.y = - sin_pitch_disk * sin_theta - cos_pitch_disk * cos_theta;

    double bMag = cos(theta - cos_pitch_disk / sin_pitch_disk * log(r / R_sun) + PHI);
    bMag *= B0_disk * R_sun / std::max(r, Rc_disk) / cos_PHI;
    B *= bMag;

    return B;
}

Vector3d FMO20Field::getHaloField(const Vector3d& pos) const {
    Vector3d B(0.);

    double rho = pos.getR(); // spherical radius

    if (rho >= rho0_halo){ 
	    double theta = pos.getTheta();
	    double phi = pos.getPhi();
	    double cos_phi = cos(phi);
	    double sin_phi = sin(phi);
	    double cos_theta = cos(theta);
	    double sin_theta = sin(theta);

        // Define using a Parker / Archimedean spiral field
	    // radial direction
	    double C1 = rho0_halo * rho0_halo / (rho * rho);
	    B.x += C1 * cos_phi * sin_theta;
	    B.y += C1 * sin_phi * sin_theta;
	    B.z += C1 * cos_theta;
	    
	    // azimuthal direction	
	    double C2 = - (rho0_halo * rho0_halo * sin_theta) / (rho * rho1_halo);
	    B.x += C2 * (-sin_phi);
	    B.y += C2 * cos_phi;

	    // magnetic field switch at z = 0
	    if (pos.z<0.) {
	    	B *= -1;
	    }

	    B *= B0_halo;
    }
    return B;
}

double FMO20Field::logisticalFunction(const double& z, const double& z0, const double& k) const {
    // to avoid mixing halo and disk 
    // -> modified logistical function acting as a smoothed Heaviside
    return 1. / (1. + exp(-k * (abs(z) - z0) )) ;

}

void FMO20Field::setB0Disk(double B) { B0_disk = B; }
void FMO20Field::setRcDisk(double Rc) { Rc_disk = Rc; }
void FMO20Field::setDDisk(double d) { d_disk = d; setParams(); }
void FMO20Field::setPitchAngle(double pitch_angle) { 
    pitch_disk = pitch_angle; 
    sin_pitch_disk = sin(pitch_angle);
    cos_pitch_disk = cos(pitch_angle);
    setParams();
}
void FMO20Field::setParams() {
	PHI = cos_pitch_disk / sin_pitch_disk * log1p(d_disk / R_sun) - M_PI / 2.;
	cos_PHI = cos(PHI + M_PI); // Normalisation to have B0 = B(Earth)
}


void FMO20Field::setB0Halo(double B) { B0_halo = B; }
void FMO20Field::setRho0Halo(double rho) { rho0_halo = rho; }
void FMO20Field::setRho1Halo(double rho) { rho1_halo = rho; }
void FMO20Field::setZ0(double z) { z0 = z; }
void FMO20Field::setKLF(double k) { k_LF = k; }

void FMO20Field::setUseDisk(bool use) { use_disk = use; }
void FMO20Field::setUseHalo(bool use) { use_halo = use; }
bool FMO20Field::isUsingDisk() { return use_disk; }
bool FMO20Field::isUsingHalo() { return use_halo; }

} // namespace crpropa
