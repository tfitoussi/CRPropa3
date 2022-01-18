#include "crpropa/magneticField/MyGMF.h"
#include "crpropa/Units.h"
#include <math.h>

namespace crpropa {

MyGMF::MyGMF(
        const double B0, // muG
        const double Rb, // kpc
        const double z0, // kpc
        const double d_, // kpc
        const double pitch, // rad
        const double rho0, // kpc
        const double rho1, // kpc
        const double lambdar, // kpc
        const double lambdaz // kpc
){
	// disk parameters
	R_b = Rb;
	d = d_;
	z_0 = z0;
	B_0 = B0;

    rho_0 = rho0;
    rho_1 = rho1;

	R_sun = 8.5 * kpc;

    lambda_gc = 1. * pc;
    lambda_r = lambdar;
    lambda_z = lambdaz;

	cos_pitch = cos(pitch);
	sin_pitch = sin(pitch);
	PHI = cos_pitch / sin_pitch * log1p(d / R_sun) - M_PI / 2;
	cos_PHI = cos(PHI);

    setUseDisk(true);
    setUseHalo(true);
}

double MyGMF::logisticalFunction(const double &x, const double &x0, const double &lambda) const {
    
    return 1. / (1. + exp(-(abs(x) - x0)/lambda)) ;
}

double MyGMF::transition(const double r, const double z) const {
    double lf_gc = logisticalFunction(r, rho_0, lambda_gc); // galactic center
    double lf_ed = logisticalFunction(r, R_b, lambda_r); // edge disk
    double lf_dh = logisticalFunction(z, z_0, lambda_z); // transistion disk - halo
    return lf_gc * (1. - lf_ed) * (1. - lf_dh);
}

Vector3d MyGMF::getField(const Vector3d& pos) const {
	double r = sqrt(pos.x * pos.x + pos.y * pos.y);  // in-plane radius
    double trans = transition(r, pos.z);

	Vector3d B(0.);

    if (use_disk) {
	    // PT11 paper has B_theta = B * cos(p) but this seems because they define azimuth clockwise, while we have anticlockwise.
	    // see Tinyakov 2002 APh 18,165: "local field points to l=90+p" so p=-5 deg gives l=85 and hence clockwise from above.
	    // so to get local B clockwise in our system, need minus (like Sun etal).
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
	    B.x = sin_pitch * cos_theta - cos_pitch * sin_theta;
	    B.y = - sin_pitch * sin_theta - cos_pitch * cos_theta;
	    B *= -1;	// flip magnetic field direction, as B_{theta} and B_{phi} refering to 180 degree rotated field

		double bMag = cos(theta - cos_pitch / sin_pitch * log(r / R_sun) + PHI);
		bMag *= B_0 * trans;
		B *= bMag;
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
	    B.x += cos_phi * sin_theta;
	    B.y += sin_phi * sin_theta;
	    B.z += cos_theta;
	    
	    // azimuthal direction	
	    double C2 = - (rho * sin_theta) / rho_1;
	    B.x += C2 * (-sin_phi);
	    B.y += C2 * cos_phi;

	    // magnetic field switch at z = 0
	    if (pos.z>0.) {
	        B *= -1;
	    }

        // constant field below rho_0, 1/rho**2 else
        double norm = rho <= rho_0 ? 1 : rho_0 * rho_0 / (rho * rho);
	    B *= B_0 * norm * (1. - trans);
    }
    return B;
}

void MyGMF::setUseDisk(bool use) { use_disk = use; }
void MyGMF::setUseHalo(bool use) { use_halo = use; }

} // namespace crpropa
