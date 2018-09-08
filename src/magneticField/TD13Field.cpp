#include "crpropa/magneticField/TD13Field.h"
#include "crpropa/Units.h"
#include "crpropa/GridTools.h"
#include "crpropa/Random.h"

#include <iostream>

#include <immintrin.h>
#include <sleef.h>

#define USE_SIMD

namespace crpropa {

std::vector<double> logspace(double start, double stop, size_t N) {

  double delta = stop - start;
  std::vector<double> values = std::vector<double>(N, 0.);
  for (int i=0; i<N; i++) {
    values[i] = pow(10, ((double) i) / ((double) (N-1)) * delta + start);
  }
  return values;
}

// code from:
// https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-float-vector-sum-on-x86
float hsum_float_sse3(__m128 v) {
    __m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
    __m128 sums = _mm_add_ps(v, shuf);
    shuf        = _mm_movehl_ps(shuf, sums); // high half -> low half
    sums        = _mm_add_ss(sums, shuf);
    return        _mm_cvtss_f32(sums);
}

  TD13Field::TD13Field(double Brms, double kmin, double kmax, double gamma, double bendoverScale, int Nm, int seed) {

    // NOTE: the use of the turbulence bend-over scale in the TD13 paper is quite confusing to
    // me. The paper states that k = l_0 * <k tilde> would be used throughout, yet
    // surely they do not mean to say that l_0 * <k tilde> should be used for the k in the
    // scalar product in eq. 2? In this implementation, I've only multiplied in the l_0
    // in the computation of the Gk, not the actual <k>s used for planar wave evaluation,
    // since this would yield obviously wrong results...

    if (kmin > kmax) {
      throw std::runtime_error("TD13Field: kmin > kmax");
    }

    if (Nm <= 1) {
      throw std::runtime_error("TD13Field: Nm <= 1. We need at least two wavemodes in order to generate the k distribution properly, and besides -- *what are you doing?!*");
    }

    if (kmin <= 0) {
      throw std::runtime_error("TD13Field: kmin <= 0");
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
    //on second thought, this is probably unnecessary since it's just a factor and will get
    //normalized out anyways.

    double Ak2_sum = 0; // sum of Ak^2 over all k
    //for this loop, the Ak array actually contains Gk*delta_k (ie non-normalized Ak^2)
    for (int i=0; i<Nm; i++) {
      double k = this->k[i] * bendoverScale;
      double Gk = pow(k, q) / pow(1 + k*k, (s+q)/2);
      Ak[i] = Gk * delta_k0 * k  *k*k; //DEBUG volume correction factor
      Ak2_sum += Ak[i];
    }
    //only in this loop are the actual Ak computed and stored
    //(this two-step process is necessary in order to normalize the values properly)
    for (int i=0; i<Nm; i++) {
      Ak[i] = sqrt(Ak[i] / Ak2_sum * 2) * Brms;
    }

    // generate direction, phase, and polarization for each wavemode
    for (int i=0; i<Nm; i++) {
      double phi = random.randUniform(-M_PI, M_PI);
      double costheta = random.randUniform(-1., 1.);
      //// DEBUG set these to zero for aligned FFT
      //phi = 0.;
      //costheta = 0.;
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
    //copy data into AVX-compatible arrays
    avx_Nm = ( (Nm + 8 - 1)/8 ) * 8; //round up to next larger multiple of 8: align is 256 = 8 * sizeof(float) bit
    //std::cout << avx_Nm <<std::endl;
    //std::cout << itotal << std::endl;
    //std::cout << (itotal*avx_Nm + 7) << std::endl;
    avx_data = std::vector<float>(itotal*avx_Nm + 7, 0.);

    //get the first 256-bit aligned element
    size_t size = avx_data.size()*sizeof(float);
    void *pointer = avx_data.data();
    align_offset = (float *) std::align(32, 32, pointer, size) - avx_data.data();

    //copy
    for (int i=0; i<Nm; i++) {
      avx_data[i + align_offset + avx_Nm*ixi0] = xi[i].x;
      avx_data[i + align_offset + avx_Nm*ixi1] = xi[i].y;
      avx_data[i + align_offset + avx_Nm*ixi2] = xi[i].z;

      avx_data[i + align_offset + avx_Nm*ikappa0] = kappa[i].x;
      avx_data[i + align_offset + avx_Nm*ikappa1] = kappa[i].y;
      avx_data[i + align_offset + avx_Nm*ikappa2] = kappa[i].z;

      avx_data[i + align_offset + avx_Nm*iAk] = Ak[i];
      avx_data[i + align_offset + avx_Nm*ik] = k[i];
      avx_data[i + align_offset + avx_Nm*ibeta] = beta[i];
    }
}

Vector3d TD13Field::getField(const Vector3d& pos) const {

#ifndef USE_SIMD
  Vector3d B(0.);
  
  for (int i=0; i<Nm; i++) {
    double z_ = pos.dot(kappa[i]);
    B += xi[i] * Ak[i] * cos(k[i] * z_ + beta[i]);
  }

  return B;

#else
  __m128 acc0 = _mm_setzero_ps();
  __m128 acc1 = _mm_setzero_ps();
  __m128 acc2 = _mm_setzero_ps();

  __m128 pos0 = _mm_set1_ps(pos.x);
  __m128 pos1 = _mm_set1_ps(pos.y);
  __m128 pos2 = _mm_set1_ps(pos.z);

  __m128 test;

  for (int i=0; i<avx_Nm; i+=4) {

    //load data from memory into AVX registers
    __m128 xi0 = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*ixi0);
    __m128 xi1 = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*ixi1);
    __m128 xi2 = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*ixi2);

    __m128 kappa0 = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*ikappa0);
    __m128 kappa1 = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*ikappa1);
    __m128 kappa2 = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*ikappa2);

    __m128 Ak = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*iAk);
    __m128 k = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*ik);
    __m128 beta = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*ibeta);

    //do the computation
    __m128 z = _mm_add_ps(_mm_mul_ps(pos0, kappa0),
			      _mm_add_ps(_mm_mul_ps(pos1, kappa1),
					    _mm_mul_ps(pos2, kappa2)
					    )
			      );

    __m128 cos_arg = _mm_add_ps(_mm_mul_ps(k, z), beta);
    __m128 mag = _mm_mul_ps(Ak, Sleef_cosf4_u10(cos_arg));

    acc0 = _mm_add_ps(_mm_mul_ps(mag, xi0), acc0);
    acc1 = _mm_add_ps(_mm_mul_ps(mag, xi1), acc1);
    acc2 = _mm_add_ps(_mm_mul_ps(mag, xi2), acc2);
  }
  
  return Vector3d(hsum_float_sse3(acc0),
                  hsum_float_sse3(acc1),
                  hsum_float_sse3(acc2)
                  );
#endif
}

} // namespace crpropa
