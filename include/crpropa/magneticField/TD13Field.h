#ifndef CRPROPA_TD13FIELD_H
#define CRPROPA_TD13FIELD_H

#include <vector>
#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/Grid.h"

namespace crpropa {

std::vector<double> logspace(double start, double stop, size_t N) {

  double delta = stop - start;
  std::vector<double> values = std::vector<double>(N, 0.);
  for (int i=0; i<N; i++) {
    values[i] = pow(10, ((double) i) / ((double) (N-1)) * delta + start);
  }
  return values;
} 

/**
 @class TD13Field
 @brief Interpolation-free turbulent magnetic field based on the TD13 paper

 blah blah blah
 */
class TD13Field: public MagneticField {
private:

  std::vector<Vector3d> xi; //TODO: this is actually psi (as defined in the paper), because I'm stupid. I don't think that's a problem, but I should probably still change it.
  std::vector<Vector3d> kappa;
  std::vector<double> phi;
  std::vector<double> costheta;
  std::vector<double> beta;

  double gamma;
  double Nm;

public:
  //DEBUG put these in public so I can read them in python
  std::vector<double> Ak;
  std::vector<double> k;
  /** Constructor
      @param root mean square field strength for generated field
      @param kmin wave number of the mode with the largest wavelength to be included in the spectrum
      @param kmax wave number of the mode with the smallest wavelength to be included in the spectrum
      @param gamma spectral index
      @param Nm number of wavemodes that will be used when computing the field. A higher value will give a more accurate representation of the turbulence, but increase the runtime for getField.
      @param seed can be used to seed the random number generator used to generate the field. This works just like in initTurbulence: a seed of 0 will lead to a randomly initialized RNG.
*/
  TD13Field(double Brms, double kmin, double kmax, double gamma, int Nm, int seed = 0);

  /**
     Theoretical runtime is O(Nm).
*/
  Vector3d getField(const Vector3d& pos) const;
};

} // namespace crpropa

#endif // CRPROPA_TD13FIELD_H
