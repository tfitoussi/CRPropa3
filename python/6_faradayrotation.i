%include "crpropa/FaradayRotation.h"

%extend crpropa::FaradayRotation{
    PyObject *computeFR_numpyArray(const double distmax, const double diststep)
    {
        double* data = $self->computeFR(distmax, diststep);
        npy_intp npix = $self->getNumberOfPixels();
        npy_intp dims[1] = {npix};
        return PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, (void*)data);
    }
};
