#include "crpropa/FaradayRotation.h"
#include "crpropa/massDistribution/NE2001.h"
#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include <stdio.h>

namespace crpropa {
using namespace std;

FaradayRotation::FaradayRotation(MagneticField* GMF_, FreeElectronDensity* Ne_, const int HPorder){
    GMF = GMF_;
    _pixelization = new crpropa::Pixelization(HPorder);
    Npix = _pixelization -> getNumberOfPixels();
    Ne = Ne_;
}

double* FaradayRotation::computeFR(const double distmax, const double diststep ) const {
    double* FRmap = new double[Npix];
    double theta, phi, Ur_x, Ur_y, Ur_z;
    Vector3d pos, B;
    double ne, Bpar;
    double D_Earth = 8.5 * kpc;

    double Npix_done = 0;

    #pragma omp parallel for private(theta, phi, Ur_x, Ur_y, Ur_z, pos, B, Bpar, ne)
    for (int ipix=0; ipix < Npix; ipix++){
        _pixelization -> pix2Direction(ipix,phi,theta);
        theta = M_PI/2. - theta; // correcting with HEALpix definition in HEALpy

        Ur_x = sin(theta) * cos(phi);
        Ur_y = sin(theta) * sin(phi);
        Ur_z = cos(theta);
        FRmap[ipix] = 0;

        for (double dist=0; dist<=distmax; dist += diststep){
            pos = Vector3d(Ur_x*dist-D_Earth, Ur_y*dist, Ur_z*dist);
            B = GMF -> getField(pos) / muG;
            Bpar = B.x * Ur_x + B.y * Ur_y + B.z * Ur_z;
            ne = Ne -> getNe(pos) /ccm;

            FRmap[ipix] = FRmap[ipix] - 0.81 * Bpar * ne * diststep/pc;
        }

        #pragma omp critical(progress)
        {
        Npix_done ++;
        //cout << "\r"<<"Computing Faraday rotation [" << Npix << " pixels] ... " << int(Npix_done/Npix*100) <<" %";
        }
    }
    //cout << endl;
    return FRmap;
}

int FaradayRotation::getNumberOfPixels() const { return Npix; }

} // namespace crpropa
