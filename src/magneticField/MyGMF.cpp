#include "crpropa/magneticField/MyGMF.h"
#include "crpropa/Units.h"

namespace crpropa {

MyGMF::MyGMF(
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
){
	// disk parameters
    setUseDisk(true);
	B0_D = B0d;   // magnetic field scale
	R_sun = 8.5 * kpc;
    R_c = Rc;
	cos_pitch = cos(pitch);
	sin_pitch = sin(pitch);
    inv_tan_pitch = cos_pitch / sin_pitch; // 1 / tan(pitch)
	PHI =  log1p(d / R_sun) * inv_tan_pitch + M_PI / 2;
	cos_PHI = cos(PHI);

    // halo parameters
    setUseHalo(true);
	B0_H = B0h;   // magnetic field scale
    rho_0 = rho0;
    rho_1 = rho1;

    // transition disk - halo
    z_0 = z0; // altitude transition disk halo
    l_dh = ldh; // transition disk halo
}

double MyGMF::logisticalFunction(const double &x, const double &x0, const double &lambda) const {
    return 1. / (1. + exp(-(x - x0)/lambda)) ;
}

double MyGMF::transition(const double z) const {
    double lf_dh = logisticalFunction(fabs(z), z_0, l_dh); // transistion disk - halo
    return 1. - lf_dh;
}

Vector3d MyGMF::getField(const Vector3d& pos) const {
	double r = sqrt(pos.x * pos.x + pos.y * pos.y);  // in-plane radius
    double rho = pos.getR(); // spherical radius
    // the following is equivalent to sin(phi) and cos(phi) which is computationally slower
    double cos_phi = pos.x / r;
    double sin_phi = pos.y / r;
    double trans = transition(pos.z);

    Vector3d Bdisk(0.);
    Vector3d Bhalo(0.);
    if (use_disk) {
		double cos_phi = pos.x / r;
		double sin_phi = pos.y / r;

	    // After some geometry calculations (on whiteboard) one finds:
	    // Bx = +cos(phi) * B_r - sin(phi) * B_{theta}
	    // By = +sin(phi) * B_r + cos(phi) * B_{theta}
	    Bdisk.x = sin_pitch * cos_phi - cos_pitch * sin_phi;
	    Bdisk.y = sin_pitch * sin_phi + cos_pitch * cos_phi;

		double angle = pos.getPhi() + log(r / R_sun)  * inv_tan_pitch - PHI;
        double B_mag = B0_D * R_sun / std::max(r, R_c) / cos_PHI;

		Bdisk *= B_mag * cos(angle) * trans;
    }

    if (use_halo) {
        if (rho > 1*kpc) {
            double theta = pos.getTheta();
            double cos_theta = cos(theta);
            double sin_theta = sin(theta);
    
            // Define using a Parker / Archimedean spiral field
    	    // radial direction
    	    Bhalo.x = cos_phi * sin_theta;
    	    Bhalo.y = sin_phi * sin_theta;
    	    Bhalo.z = cos_theta;
    	    
    	    // azimuthal direction	
    	    double C2 = + rho / rho_1 * sin_theta;
    	    Bhalo.x += C2 * (-sin_phi);
    	    Bhalo.y += C2 * cos_phi;
    
    	    // magnetic field switch at z = 0
    	    if (pos.z<0.) { Bhalo *= -1; }
    
    	    Bhalo *= B0_H * rho_0*rho_0 / (rho*rho) * (1. - trans);
        }
    }

    return Bdisk + Bhalo;
}

void MyGMF::setUseDisk(bool use) { use_disk = use; }
void MyGMF::setUseHalo(bool use) { use_halo = use; }

} // namespace crpropa
