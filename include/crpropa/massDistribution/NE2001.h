#ifndef CRPROPA_NE2001
#define CRPROPA_NE2001

#include "crpropa/massDistribution/FreeElectronDensity.h"
#include "crpropa/Vector3.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"
#include <cmath>
#include <string>

namespace crpropa {
using namespace std;

/* @class   NE2001
 * @brief   NE2001 free electron density model
 *
 * Implementation of the free electron density of Cordes and Lazzio 2002,
 * NE2001. I. A NEW MODEL FOR THE GALACTIC DISTRIBUTION OF FREE ELECTRONS AND ITS FLUCTUATIONS
 * See https://arxiv.org/abs/astro-ph/0207156v3
 *
 * Model corrected with  Schnitzeler 2012,
 * Modelling the Galactic distribution of free electrons
 * doi: 10.1111/j.1365-2966.2012.21869.x
*/

class NE2001 : public FreeElectronDensity {
private:
    // Parameters Thick Disk
    double n1;
    double H1;
    double A1;

    // Parameters Thin Disk
    double n2;
    double H2;
    double A2;
    
    // Parameters Spiral Arms
    int Narms;
    double na;
    double n_arm[5];
    double w_arm[5];
    double h_arm[5];
    double alpha_arm[5];
    double R_min_arm[5];
    double theta_min_arm[5];
    double extent_arm[5];
    double Ha;
    double Aa;
    double wa;
    double theta_step_arm;
    // better performance with static array than std::vector
    double radius_arm[5][720];  // 720 = 360 deg / 0.5 deg max 
    int n_theta_arm[5]; // number of step in the array
    // to deform arm 2 and 4
    double theta_arm2_lim[4]; 
    double theta_arm2[3];
    double theta_arm4_lim[2]; 
    double theta_arm4; 
    double theta2a, theta2b; 
    double theta3a, theta3b; 
    
    // Parameters Galactic Center (GC)
    double n_gc0;
    double R_gc;
    double H_gc;
    double x_gc, y_gc, z_gc;
    
    // Parameters Voids
    int N_voids;
    double ne_void[17];
    double x_void[17];
    double y_void[17];
    double z_void[17];
    double aa_void[17];
    double bb_void[17];
    double cc_void[17];
    double c2_void[17];
    double s2_void[17];
    double c1_void[17];
    double s1_void[17];
    double c1c2_void[17];
    double s1c2_void[17];
    double c1s2_void[17];
    double s1s2_void[17];

    // Parameters Clumps
    int N_clumps;
    double ne_clumps[175];
    double x_clumps[175];
    double y_clumps[175];
    double z_clumps[175];
    double R_clumps[175];
    
    // Parameters Local Interstellar Medium (LISM)
    double apdr, bpdr, cpdr, dpdr, xldr, yldr, zldr, neldr0;
    double apsb, bpsb, cpsb, dpsb, xlsb, ylsb, zlsb, nelsb0;
    double alhb, blhb, clhb, yzslope, xlhb, ylhb, zlhb, nelhb0;
    double xlpI, ylpI, zlpI, rlpI, alpI, nelpI, dnelpI;
    
    // Common parameters
    double Rsun;

public:
    /** Construtor
     */
    NE2001();

    /**@brief   get total density of free electrons at a given position
     * @param   pos : crpropa::Vector3d position
     * @return  total free electron density in cm^-3
     */
    double getNe(const Vector3d& pos) const;

    /**@brief   get density of free electrons in the disk at a given position 
     * @param   pos : crpropa::Vector3d position
     * @return  free electron density in cm^-3 in the disk (thin disk, thick disk and arms)
     */
    double getNeDisk(const Vector3d& pos) const;
    double getNeThickDisk(const Vector3d& pos) const;
    double getNeThinDisk(const Vector3d& pos) const;
    double getNeSpiralArms(const Vector3d& pos) const;

    /**@brief   get density of free electrons near the Galactic center
     * @param   pos : crpropa::Vector3d position
     * @return  free electron density in cm^-3 
     */
    double getNeGalacticCenter(const Vector3d& pos) const;

    /**@brief   get density of free electrons if in the Clumps
     * @param   pos : crpropa::Vector3d position
     * @return  free electron density in cm^-3 
     */
    double getNeClumps(const Vector3d& pos) const;

    /**@brief   get density of free electrons if in a void 
     * @param   pos : crpropa::Vector3d position
     * @return  free electron density in cm^-3 
     */
    double getNeVoids(const Vector3d& pos) const;

    /**@brief   get density of free electrons if in the local interstellar medium
     * @param   pos : crpropa::Vector3d position
     * @return  free electron density in cm^-3 
     */
    double getNeLISM(const Vector3d& pos) const;

    /**@brief   Compute the square of the hyperbolic secant of x
     * @param   x
     * @return  ( 2 / (exp(x) + exp(-x))**2
     */
    double sech2(const double x) const;

    /**@brief   Compute the minimal distance to a given arm in the cartesian plan (x, y)
     * @param x, y      cartesian coordinates
     * @param j_arm     number of the arm 
     * @return s_min    minimal distance to the arm j
     */
    double SjX(const double x, const double y, const int j_arm) const;

    /**@brief   Generate the tables containing the cylindrical coordinates (theta, r) of each arms
     */
    void generateArms();

    /**@brief   Load datafile containing coordinates and free electrons density of the voids
     */
    void loadNeVoid();

    /**@brief   Load datafile containing coordinates and free electrons density of the local Clumps
     */
    void loadNeClumps();
};

} // CRPROPA NAMESPACE

#endif // CRPROPA_NE2001
