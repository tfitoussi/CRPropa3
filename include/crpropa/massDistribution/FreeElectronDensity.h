#ifndef CRPROPA_FED
#define CRPROPA_FED

#include "crpropa/Vector3.h"

namespace crpropa {

/* @class   FreeElectronDensity
 * @brief   Generic class for free electron density models
 *
 * Do not use directly ! Use NE2201 or YMW17 instead
*/

class FreeElectronDensity {
public:
    virtual double getNe(const Vector3d& pos) const { return 0;};
};

} // CRPROPA NAMESPACE

#endif // CRPROPA_NE2001
