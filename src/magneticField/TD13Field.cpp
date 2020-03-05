#include "crpropa/magneticField/TD13Field.h"
#include "crpropa/Units.h"
#include "crpropa/GridTools.h"
#include "crpropa/Random.h"
#include "kiss/logger.h"

#include <iostream>

#if defined(CRPROPA_HAVE_SLEEF) && defined(__SSE__) && defined(__SSE2__) && defined(__SSE3__)
#define FAST_TD13
#include <immintrin.h>
#include <sleef.h>
#endif

namespace crpropa {
#ifdef FAST_TD13
// code from:
// https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-float-vector-sum-on-x86
float hsum_float_sse3(__m128 v) {
    __m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
    __m128 sums = _mm_add_ps(v, shuf);
    shuf        = _mm_movehl_ps(shuf, sums); // high half -> low half
    sums        = _mm_add_ss(sums, shuf);
    return        _mm_cvtss_f32(sums);
}

  // TODO: Figure out why ld can't find Sleef_x86CpuID
  //  (Alternatively, decide that this isn't worth it)
  //bool check_sse() {
  //  int32_t result[4];
  //  Sleef_x86CpuID(result, 1, 0);
  //  return (result[3] & (1 << 25)) && (result[3] & (1 << 26)) && (result[2] & (1 << 0))
  //    && (result[2] & (1 << 28)) //DEBUG check for avx, to test if this is working.
	//    ;
  //}
#endif // defined(FAST_TD13)

TD13Field::TD13Field(double Brms, double Lmin, double Lmax, double s, double q, int Nm, int seed, const bool powerlaw) {

#ifdef FAST_TD13
    KISS_LOG_INFO << "TD13Field: Using SIMD TD13 implementation" << std::endl;

    // In principle, we could dynamically dispatch to the non-SIMD version in
    // this case. However, this complicates the code, incurs runtime overhead,
    // and is unlikely to happen since SSE3 is quite well supported.
    // TODO: this is currently uncommented b/c sleef seems to fail to provide
    // the cpuid function
    //if (!check_sse()) {
    //  throw std::runtime_error("TD13Field: This code was compiled with SIMD support (SSE1-3), but it is running on a CPU that does not support these instructions. Please set USE_SIMD to OFF in CMake and recompile CRPropa.");
    //}
#endif

    if (Nm <= 1) {
        throw std::runtime_error("TD13Field: Nm <= 1. We need at least two wavemodes in order to generate the k distribution properly, and besides -- *what are you doing?!*");
    }

    // initialize everything
    this->Brms = Brms;
    this->Nm = Nm;
    this->Lmax = Lmax;
    this->Lmin = Lmin;
    this->q = q;
    this->s = s;

    Random random;
    if (seed != 0) random.seed(seed);

    double kmax = 2*M_PI / Lmin;
    double kmin = 2*M_PI / Lmax;

    xi = std::vector<Vector3d>(Nm, Vector3d(0.));
    kappa = std::vector<Vector3d>(Nm, Vector3d(0.));
    phi = std::vector<double>(Nm, 0.);
    costheta = std::vector<double>(Nm, 0.);
    beta = std::vector<double>(Nm, 0.);
    Ak = std::vector<double>(Nm, 0.);
    k = std::vector<double>(Nm, 0.);

    double delta = log10(kmax/kmin);
    for (int i=0; i<Nm; i++) {
        k[i] = pow(10, log10(kmin) + ((double) i) / ((double) (Nm-1)) * delta);
    }

    // compute Ak
    double delta_k0 = (k[1] - k[0]) / k[1]; // multiply this by k[i] to get delta_k[i]
    //on second thought, this is probably unnecessary since it's just a factor and will get
    //normalized out anyways.

    double Ak2_sum = 0; // sum of Ak^2 over all k
    //for this loop, the Ak array actually contains Gk*delta_k (ie non-normalized Ak^2)
    for (int i=0; i<Nm; i++) {
        double k = this->k[i];
        double Gk;
        if (powerlaw) {
            Gk = pow(k, -s);
        } else {
            Gk = pow(k,q) / pow( 1+k*k, (s+q)/2. );
        }
        Ak[i] = Gk * delta_k0 * k;
        Ak2_sum += Ak[i];

        // phi, costheta, and sintheta are for drawing vectors with
        // uniform distribution on the unit sphere.
        // This is similar to Random::randVector(): their t is our phi,
        // z is costheta, and r is sintheta. Our kappa is equivalent to
        // the return value of randVector(); however, TD13 then reuse
        // these values to generate a random vector perpendicular to kappa.
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

    // To normalize the spectrum to Brms, we basically need to compute
    // the area under the spectral distribution. Now, if the <delta_k>s
    // are all correct, this integral is actually approximated
    // numerically by Ak2_sum. That is, at least, so long as we
    // consider the integral from kmin to kmax.
    //
    // However, from my current point of view, kmin seems to have a
    // lot more physical relevance than kmax. The correlation length
    // is clearly a physical property of the field, and it's formula
    // depends on both kmin and kmax. However, when kmin and kmax are
    // very far apart (maybe two orders of magnitude), the correlation
    // length effectively stops depending on kmax and becomes nearly
    // linear in kmin.
    //
    // So what meaning, then, does kmax have? This is actually an
    // important question, because the normalization factor for all
    // other wavemodes depends on it. Simply adding more wavemodes at
    // the top will scale all other wavemodes down, to conserve
    // Brms. The effect of this will not be particularly large,
    // considering that the higher wavemodes only contribute very
    // little due to the power law distribution.
    //
    // This leads to the approach of normalizing the wavemodes as if
    // the spectrum went from kmin to infinity. This removes kmax
    // completely from the equation, but the additional contribution
    // stays negligible.
    //
    // ## Implementation
    //
    // So, instead of dividing by Ak2_sum (which is roughly equivalent
    // to the integral from kmin to kmax), we want to normalize with
    // the integral from kmin to infinity. We could simply divide by
    // this integral instead, but then the correctness of the
    // resulting value would depend on the correctness of the
    // <delta_k>s, and I don't really want to rely on that. The more
    // elegant solution seems to be (to me), to divide out the Ak2_sum
    // first, nullifying any problems with the <delta_k>s, then
    // multiplying the value from kmin to kmax back in, and finally
    // normalizing to the value of the integral from kmin to infinity.

    // integral from kmin to kmax
    double theoretical_Ak2_sum = 1/(s-1) * (pow(kmin, -s+1) - pow(kmax, -s+1));
    // integral from kmin to infinity
    double normalization_factor = 1/(s-1) * pow(kmin, -s+1);

    Ak2_sum = Ak2_sum * theoretical_Ak2_sum / normalization_factor;
    
    // only in this loop are the actual Ak computed and stored
    // (this two-step process is necessary in order to normalize the values properly)
    for (int i=0; i<Nm; i++) {
        Ak[i] = sqrt(2 * Ak[i] / Ak2_sum ) * Brms;
    }


#ifdef FAST_TD13
    // copy data into AVX-compatible arrays
    //
    // What is going on here:
    //
    // SIMD load instructions require data to be aligned in memory to
    // the SIMD register bit size (128 bit for SSE, 256 bit for
    // AVX). Specifying memory alignment in C++ felt somewhat icky to
    // me after doing initial research, so here we're using the simple
    // solution: Allocating a vector that is a little bit longer than
    // required, determining the first aligned index, and then only
    // using the array starting at that index. (Here called align_offset.)
    //
    // To cut down on the complexity of storing one such offset for
    // each of the attribute arrays, we're using only a single array to
    // store all of the wavemode attributes, with offset indices to
    // each. This way, only a single align_offset has to be stored.
    //
    // This code was originally written for AVX, hence the
    // naming. I've also kept the larger 256-bit alignment required
    // for AVX, even though SSE only needs 128-bit alignment, to make
    // it simpler to go back to AVX. (And because there's really no
    // disadvantage, except for allocating a few bytes more.)
    //
    // The second thing to consider is that SSE reads data in chunks
    // of four floats. So what happens when the number of modes is not
    // divisible by four? For this, we need to round up the number of
    // modes to the next multiple of four. (Again, this code rounds up
    // to eight, due to AVX compatibility. TODO: maybe remove this,
    // since it actually costs runtime?)
    //
    // Since the array is initialized to zero, the "extra modes" will
    // have all their attributes be zero. The final step in the
    // getField() loop is multiplying by (Ak*xi), which is zero for
    // the extra modes, so they won't influence the result.
    
    avx_Nm = ( (Nm + 8 - 1)/8 ) * 8; //round up to next larger multiple of 8: align is 256 = 8 * sizeof(float) bit
    avx_data = std::vector<float>(itotal*avx_Nm + 7, 0.);

    //get the first 256-bit aligned element
    size_t size = avx_data.size()*sizeof(float);
    void *pointer = avx_data.data();
    align_offset = (float *) std::align(32, 32, pointer, size) - avx_data.data();

    //copy
    for (int i=0; i<Nm; i++) {
        avx_data[i + align_offset + avx_Nm*iAxi0] = Ak[i] * xi[i].x;
        avx_data[i + align_offset + avx_Nm*iAxi1] = Ak[i] * xi[i].y;
        avx_data[i + align_offset + avx_Nm*iAxi2] = Ak[i] * xi[i].z;

        avx_data[i + align_offset + avx_Nm*ikkappa0] = k[i] * kappa[i].x;
        avx_data[i + align_offset + avx_Nm*ikkappa1] = k[i] * kappa[i].y;
        avx_data[i + align_offset + avx_Nm*ikkappa2] = k[i] * kappa[i].z;

        avx_data[i + align_offset + avx_Nm*ibeta] = beta[i];
    }
#endif
}

Vector3d TD13Field::getField(const Vector3d& pos) const {
#ifndef FAST_TD13
    Vector3d B(0.);
    for (int i=0; i<Nm; i++) {
        double z_ = pos.dot(kappa[i]);
        B += xi[i] * Ak[i] * cos(k[i] * z_ + beta[i]);
    }
    return B;

#else // CRPROPA_USE_SIMD
    // Initialize accumulators
    //
    // There is one accumulator per component of the result vector.
    // Note that each accumulator contains four numbers. At the end of
    // the loop, each of these number will contain the sum of every
    // fourth wavemodes, starting at a different offset. In the end, all
    // of the accumulator's numbers are added together (using
    // hsum_float_sse3), resulting in the total sum.

    __m128 acc0 = _mm_setzero_ps();
    __m128 acc1 = _mm_setzero_ps();
    __m128 acc2 = _mm_setzero_ps();

    // broadcast position into SSE registers
    __m128 pos0 = _mm_set1_ps(pos.x);
    __m128 pos1 = _mm_set1_ps(pos.y);
    __m128 pos2 = _mm_set1_ps(pos.z);

    for (int i=0; i<avx_Nm; i+=4) {
        // load data from memory into AVX registers
        __m128 Axi0 = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*iAxi0);
        __m128 Axi1 = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*iAxi1);
        __m128 Axi2 = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*iAxi2);

        __m128 kkappa0 = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*ikkappa0);
        __m128 kkappa1 = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*ikkappa1);
        __m128 kkappa2 = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*ikkappa2);

        __m128 beta = _mm_load_ps(avx_data.data() + i + align_offset + avx_Nm*ibeta);


        // Do the computation

        // this is the scalar product between k*kappa and pos
        __m128 z = _mm_add_ps(_mm_mul_ps(pos0, kkappa0),
        		      _mm_add_ps(_mm_mul_ps(pos1, kkappa1),
        				    _mm_mul_ps(pos2, kkappa2)
        				    )
        		      );

        // here, the phase is added on. this is the argument of the cosine.
        __m128 cos_arg = _mm_add_ps(z, beta);
        // the result of the cosine
        __m128 mag = Sleef_cosf4_u35(cos_arg);

        // Finally, Ak*xi is multiplied on. Since this is a vector, the
        // multiplication needs to be done for each of the three
        // components, so it happens separately.
        acc0 = _mm_add_ps(_mm_mul_ps(mag, Axi0), acc0);
        acc1 = _mm_add_ps(_mm_mul_ps(mag, Axi1), acc1);
        acc2 = _mm_add_ps(_mm_mul_ps(mag, Axi2), acc2);
    }
    
    return Vector3d(hsum_float_sse3(acc0),
                    hsum_float_sse3(acc1),
                    hsum_float_sse3(acc2)
                    );
#endif // FAST_TD13
}

double TD13Field::getLc() const {
    // According to Harari et Al JHEP03(2002)045
    double Lc;
    Lc = Lmax/2.;
    Lc*= (s-1.)/s;
    Lc*= 1 - pow(Lmin/Lmax, s);
    Lc/= 1 - pow(Lmin/Lmax, s-1);
    return Lc;
}

} // namespace crpropa
