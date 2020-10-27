#ifndef CRPROPA_YMW17
#define CRPROPA_YMW17

#include "crpropa/massDistribution/FreeElectronDensity.h"
#include "crpropa/Vector3.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"
#include <cmath>
#include <string>

namespace crpropa {
using namespace std;

/* @class   YMW17
 * @brief   YMW17 free electron density model
 *
 * Implementation of the free electron density of Cordes and Lazzio 2002,
 * YMW17. I. A NEW MODEL FOR THE GALACTIC DISTRIBUTION OF FREE ELECTRONS AND ITS FLUCTUATIONS
 * See https://arxiv.org/abs/astro-ph/0207156v3
 *
 * Model corrected with  Schnitzeler 2012,
 * Modelling the Galactic distribution of free electrons
 * doi: 10.1111/j.1365-2966.2012.21869.x
*/

class YMW17 : public FreeElectronDensity {
private:
    // Common parameters
    double Rsun;
    double gamma_w, phi_w, R_w; // Galactic warp

    // Parameters Thick Disk
    double Ad, Bd;
    double n1, H1;

    // Parameters Fermi Bubbles
    double J_FB, z_FB, a_FB, b_FB;

    // Parameters Thin Disk
    double n2, K2;
    double A2, B2;
    
    // Parameters Spiral Arms
    int Narms;
    double n0_arm[5];
    double w_arm[5];
    double R0_arm[5];
    double phi0_arm[5];
    double psi0_arm[5];
    double Ka, Aa;
    double nCN, phiCN, delta_phiCN;
    double nSG, phiSG, delta_phiSG;
    double phi_step_arm;
    int N0_phi;
    // better performance with static array than std::vector
    double radius_arm[5][720];  // 720 = 2 * 360 deg / 1 deg max 
    double n_arm[5][720];  // 720 = 360 deg / 1 deg max 
 
    // Parameters Galactic Center (GC)
    double x_GC, y_GC, z_GC;
    double n_GC;
    double A_GC, H_GC;

    // Parameters Gum nebula
    double x_GN, y_GN, z_GN;
    double A_GN, A_GN_2;
    double K_GN, KGN_AGN_2;
    double n_GN;
    double W_GN, W_GN_2;

    // Parameters Local bubble
    double R_LB, J_LB;
    double n_LB1, l_LB1, delta_l_LB1, W_LB1, H_LB1;
    double n_LB2, l_LB2, delta_l_LB2, W_LB2, H_LB2;

    // Parameters Loop I
    double n_LI, R_LI, W_LI;
    double l_LI, cos_l_LI, sin_l_LI, delta_l_LI;
    double x_LI, y_LI, z_LI;


    
public:
    /** Construtor
     */
    YMW17();

    /**@brief   get total density of free electrons at a given position
     * @param   pos : crpropa::Vector3d position
     * @return  total free electron density in cm^-3
     */
    double getNe(const Vector3d& pos) const;

    /**@brief   get density of free electrons in the disk at a given position 
     * @param   pos : crpropa::Vector3d position
     * @return  free electron density in cm^-3 in the disk (thin disk, thick disk and arms)
     */
    double getNeThickDisk(const Vector3d& pos) const;
    double getNeThinDisk(const Vector3d& pos) const;
    double getNeSpiralArms(const Vector3d& pos) const;

    /**@brief   get density of free electrons near the Galactic center
     * @param   pos : crpropa::Vector3d position
     * @return  free electron density in cm^-3 
     */
    double getNeGalacticCenter(const Vector3d& pos) const;

    /**@brief   get density of free electrons if in the Gum nebula
     * @param   pos : crpropa::Vector3d position
     * @return  free electron density in cm^-3 
     */
    double getNeGumNebula(const Vector3d& pos) const;

    /**@brief   get density of free electrons if in the local bubbles (1 and 2)
     * @param   pos : crpropa::Vector3d position
     * @return  free electron density in cm^-3 
     */
    double getNeLocalBubbles(const Vector3d& pos) const;

    /**@brief   get density of free electrons if in the Loop I
     * @param   pos : crpropa::Vector3d position
     * @return  free electron density in cm^-3 
     */
    double getNeLoopI(const Vector3d& pos) const;

    double sech2(const double x) const;
    double zw(const double R, const double phi) const;
    double gd(const double R) const;
    double rLB(const double x, const double y, const double z) const;
    double Hfunction(const double R) const;

    /**@brief   Generate the tables containing the cylindrical coordinates (phi, r) of each arms
     */
    void generateArms();


};

} // CRPROPA NAMESPACE

#endif // CRPROPA_YMW17
