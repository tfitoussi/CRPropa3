#include "crpropa/magneticField/MyGMF.h"
#include "crpropa/Units.h"

namespace crpropa {

MyGMF::MyGMF(
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
){
	// disk parameters
    setUseDisk(true);
	B_0d = B0d;   // magnetic field scale
	cos_pitch = cos(pitch);
	sin_pitch = sin(pitch);
    tan_pitch = sin_pitch / cos_pitch;
	R_sun = 8.5 * kpc;
    Phi_0 = phi0;

    // halo parameters
    setUseHalo(true);
	B_0h = B0h;   // magnetic field scale
    rho_0 = rho0;
    rho_1 = rho1;

    // transition disk - halo
    R_gc = Rgc; // radius central region
    l_gc = lgc; // transition central region / disk
    R_b = Rb; // scale exponential decay disk at Earth
    z_0 = z0; // altitude transition disk halo
    l_dh = ldh; // transition disk halo
}

double MyGMF::logisticalFunction(const double &x, const double &x0, const double &lambda) const {
    return 1. / (1. + exp(-(x - x0)/lambda)) ;
}

double MyGMF::transition(const double r, const double z) const {
    double lf_gc = logisticalFunction(r, R_gc, l_gc); // galactic center
    double trans_ed =  exp(-(r - R_sun)/R_b); // edge disk
    double lf_dh = logisticalFunction(fabs(z), z_0, l_dh); // transistion disk - halo
    return lf_gc * trans_ed * (1. - lf_dh);
}

Vector3d MyGMF::getField(const Vector3d& pos) const {
	double r = sqrt(pos.x * pos.x + pos.y * pos.y);  // in-plane radius
    double trans = transition(r, pos.z);

    Vector3d Bdisk(0.);
    Vector3d Bhalo(0.);

    if (use_disk) {
	    // the following is equivalent to sin(phi) and cos(phi) which is computationally slower
	    double cos_phi = pos.x / r;
	    double sin_phi = pos.y / r;

	    // After some geometry calculations (on whiteboard) one finds:
	    // Bx = +cos(theta) * B_r - sin(theta) * B_{theta}
	    // By = -sin(theta) * B_r - cos(theta) * B_{theta}
	    Bdisk.x = sin_pitch * cos_phi - cos_pitch * sin_phi;
	    Bdisk.y = sin_pitch * sin_phi + cos_pitch * cos_phi;

		double angle = pos.getPhi() - log(r / R_sun) / tan_pitch + Phi_0;
		Bdisk *= B_0d * trans * cos(angle);
    } 

    if (use_halo) {
        double rho = pos.getR(); // spherical radius
        double theta = pos.getTheta();
        double phi = pos.getPhi();
        double cos_phi = cos(phi);
        double sin_phi = sin(phi);
        double cos_theta = cos(theta);
        double sin_theta = sin(theta);

        // Define using a Parker / Archimedean spiral field
	    // radial direction
	    Bhalo.x = cos_phi * sin_theta;
	    Bhalo.y = sin_phi * sin_theta;
	    Bhalo.z = cos_theta;
	    
	    // azimuthal direction	
	    double C2 = - (rho * sin_theta) / rho_1;
	    Bhalo.x += C2 * (-sin_phi);
	    Bhalo.y += C2 * cos_phi;

	    // magnetic field switch at z = 0
	    if (pos.z>0.) {
	        Bhalo *= -1;
	    }

        // constant field below rho_0, 1/rho**2 else
        double scale = rho <= R_gc ?  R_gc*R_gc : rho*rho;
	    Bhalo *= B_0h * rho_0 * rho_0 / scale * (1. - trans);
    }
    return Bdisk + Bhalo;
}

void MyGMF::setUseDisk(bool use) { use_disk = use; }
void MyGMF::setUseHalo(bool use) { use_halo = use; }

} // namespace crpropa
