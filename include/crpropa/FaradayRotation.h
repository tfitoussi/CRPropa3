#ifndef CRPROPA_FARADAYROTATION_H
#define CRPROPA_FARADAYROTATION_H
#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/magneticLens/Pixelization.h"
#include "crpropa/massDistribution/FreeElectronDensity.h"


namespace crpropa {
/**
 * \addtogroup Core
 * @{
 */

/**
 @class FaradayRotation
 @brief Compute the HEALpix Faraday rotation map with a given CRPropa Galactic magnetic field
 */
class FaradayRotation: public Referenced {
private:
    MagneticField* GMF;
    Pixelization* _pixelization;
    int Npix;
    FreeElectronDensity* Ne;

public:
    FaradayRotation(MagneticField* GMF_, FreeElectronDensity* Ne_, const int HPorder);
    double* computeFR(const double distmax, const double diststep) const;

    int getNumberOfPixels() const;

};

} // namespace crpropa

#endif // CRPROPA_FARADAYROTATION_H
