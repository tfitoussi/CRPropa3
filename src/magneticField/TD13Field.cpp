#include "crpropa/magneticField/TD13Field.h"
#include "crpropa/Units.h"
#include "crpropa/GridTools.h"
#include "crpropa/Random.h"

#include <iostream>

namespace crpropa {

  TD13Field::TD13Field(double Brms, double kmin, double kmax, double gamma, int Nm, int seed) {

    if (kmin > kmax) {
      throw std::runtime_error("TD13Field: kmin > kmax");
    }

    if (Nm <= 1) {
      throw std::runtime_error("TD13Field: Nm <= 1. We need at least two wavemodes in order to generate the k distribution properly, and besides -- *what are you doing?!*");
    }

    if (kmin < std::numeric_limits<double>::epsilon()) {
      throw std::runtime_error("TD13Field: kmin ~<= 0");
    } 

    Random random;
    if (seed != 0) { // copied from initTurbulence
      random.seed(seed);
    }

    // initialize everything
    this->gamma = gamma;
    this->Nm = Nm;

    xi = std::vector<Vector3d>(Nm, Vector3d(0.));
    kappa = std::vector<Vector3d>(Nm, Vector3d(0.));
    phi = std::vector<double>(Nm, 0.);
    costheta = std::vector<double>(Nm, 0.);
    beta = std::vector<double>(Nm, 0.);
    Ak = std::vector<double>(Nm, 0.);

    k = logspace(log10(kmin), log10(kmax), Nm);

    // compute Ak
    double q = 0; // TODO: what is q
    double s = gamma;
    double delta_k0 = (k[1] - k[0]) / k[1]; // multiply this by k[i] to get delta_k[i]

    double Ak2_sum = 0; // sum of Ak^2 over all k
    //for this loop, the Ak array actually contains Gk*delta_k (ie non-normalized Ak^2)
    for (int i=0; i<Nm; i++) {
      double Gk = pow(k[i], q) / pow(1 + k[i]*k[i], (s+q)/2);
      Ak[i] = Gk * delta_k0 * k[i];
      Ak2_sum += Ak[i];
    }
    //only in this loop are the actual Ak computed and stored
    //(this two-step process is necessary in order to normalize the values properly)
    for (int i=0; i<Nm; i++) {
      Ak[i] = sqrt(Ak[i] / Ak2_sum * 2) * Brms;
    }

    // generate direction, phase, and polarizarion for each wavemode
    for (int i=0; i<Nm; i++) {
      double phi = random.randUniform(-M_PI, M_PI);
      double costheta = random.randUniform(-1., 1.);
      double sintheta = sqrt(1 - costheta*costheta);

      double alpha = random.randUniform(0, 2*M_PI);
      double beta = random.randUniform(0, 2*M_PI);

      Vector3d kappa = Vector3d ( sintheta * cos(phi), sintheta*sin(phi), costheta );
      Vector3d xi = Vector3d ( costheta*cos(phi)*cos(alpha) + sin(phi)*sin(alpha),
		     costheta*sin(phi)*cos(alpha) - cos(phi)*sin(alpha),
		     -sintheta*cos(alpha) );

      this->xi[i] = xi;
      this->kappa[i] = kappa;
      this->phi[i] = phi;
      this->costheta[i] = costheta;
      this->beta[i] = beta;
    }
}

Vector3d TD13Field::getField(const Vector3d& pos) const {
  Vector3d B(0.);

  for (int i=0; i<Nm; i++) {
    double z_ = pos.dot(kappa[i]);
    B += xi[i] * Ak[i] * cos(k[i] * z_ + beta[i]);
  }
  return B;
}

} // namespace crpropa
