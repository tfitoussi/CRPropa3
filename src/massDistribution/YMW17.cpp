#include "crpropa/massDistribution/YMW17.h"
#include <algorithm>
#include <fstream>

namespace crpropa {
using namespace std;

YMW17::YMW17(){
    // All densities are initially in /ccm, we multiply them by ccm to have them without dimension
    // So you can change the output by dividing it by ccm or m^3 for example

    // Common parameters
    Rsun = 8.3 * kpc;
    gamma_w = 0.14;
    phi_w = 0. * deg;
    R_w = 8.4 * kpc;
    
    // Parameters Thick Disk
    Ad = 2.5 * kpc;
    Bd = 15 * kpc;
    n1 = 0.01132 * ccm;
    H1 = 1.673 * kpc;

    // Parameters Thin Disk
    A2 = 1.2 * kpc;
    B2 = 4 * kpc;
    n2 = 0.404 * ccm;
    K2 = 1.54;
    
    // Parameters Fermi Bubbles
    J_FB = 1.;
    z_FB = 0.5 * Rsun * tan(50*deg);
    a_FB = 0.5 * Rsun * tan(50*deg);
    b_FB = Rsun * tan(20*deg);
    
    // Parameters Spiral Arms
    Narms = 5;
    phi_step_arm = 1 * deg;
    N0_phi = 720;
    n0_arm[0] = 0.135 *ccm;
    n0_arm[1] = 0.129 *ccm;
    n0_arm[2] = 0.103 *ccm;
    n0_arm[3] = 0.116 *ccm;
    n0_arm[4] = 0.0057 *ccm;
    w_arm[0] = 300 * pc;
    w_arm[1] = 500 * pc;
    w_arm[2] = 300 * pc;
    w_arm[3] = 500 * pc;
    w_arm[4] = 300 * pc;
    R0_arm[0] = 3.35 * kpc;
    R0_arm[1] = 3.71 * kpc;
    R0_arm[2] = 3.56 * kpc;
    R0_arm[3] = 3.67 * kpc;
    R0_arm[4] = 8.21 * kpc;
    phi0_arm[0] =  44.4 * deg - M_PI; // need to rotate the arms to work properly
    phi0_arm[1] = 120.0 * deg - M_PI; // need to rotate the arms to work properly
    phi0_arm[2] = 218.6 * deg - M_PI; // need to rotate the arms to work properly
    phi0_arm[3] = 330.3 * deg - M_PI; // need to rotate the arms to work properly
    phi0_arm[4] =  55.1 * deg - M_PI; // need to rotate the arms to work properly
    psi0_arm[0] = 11.43 * deg;
    psi0_arm[1] =  9.84 * deg;
    psi0_arm[2] = 10.38 * deg;
    psi0_arm[3] = 10.54 * deg;
    psi0_arm[4] =  2.77 * deg;

    Aa = 11.68 * kpc;
    Ka = 5.01;
    nCN = 2.4;
    phiCN = 109 * deg + M_PI; // need to rotate the arms to work properly
    delta_phiCN = 8.2 * deg;
    nSG = 0.626;
    phiSG = 75.8 *deg + M_PI; // need to rotate the arms to work properly
    delta_phiSG = 20 * deg;

    generateArms();

    // Parameters Galactic Center (GC)
    x_GC = 50 * pc;
    y_GC = 0 * pc;
    z_GC = -7 * pc;
    n_GC = 6.2 * ccm;
    A_GC = 160 * pc;
    H_GC = 35 * pc;

    // Parameters Local bubble
    J_LB = 0.48;
    R_LB = 110 * pc; 
    n_LB1 = 1.094 * ccm;
    J_LB = 0.48;
    l_LB1 = 195.4 * deg - M_PI; // need to rotate the arms to work properly; 
    delta_l_LB1 = 28.4 * deg; 
    W_LB1 = 14.2 * pc; 
    H_LB1 = 112.9 * pc;
    n_LB2 = 2.33 * ccm;
    l_LB2 = 278.2 * deg - M_PI; // need to rotate the arms to work properly; 
    delta_l_LB2 = 14.7 * deg; 
    W_LB2 = 15.6 * pc; 
    H_LB2 = 43.6 * pc;

    // Parameters Gum nebula
    x_GN = -446 * pc;
    y_GN = Rsun + 47*pc;
    z_GN = -31 * pc;
    n_GN = 1.84 * ccm;
    W_GN = 15.1 * pc;
    A_GN = 125.8 * pc;
    K_GN = 1.4;
    A_GN_2 = A_GN * A_GN;
    W_GN_2 = W_GN * W_GN;

    // Parameters Loop I
    n_LI = 1.907 * ccm;
    R_LI = 80 * pc; 
    W_LI = 15 * pc;
    delta_l_LI = 30 * deg;
    l_LI = 40 *deg;
    cos_l_LI = cos(l_LI); 
    sin_l_LI = sin(l_LI); 
    x_LI = 10 * pc; 
    y_LI = 8106 * pc; 
    z_LI = 10 * pc;

}

double YMW17::getNe(const Vector3d& pos) const  {
    double x =  pos.getY(); // reverse axis
    double y = -pos.getX();
    double z =  pos.getZ();

    double Ne_Local_Bubbles = getNeLocalBubbles(pos);
    double Ne_thick_disk = getNeThickDisk(pos);
    double Ne_thin_disk = getNeThinDisk(pos);
    double Ne_spiral_arms = getNeSpiralArms(pos);
    double Ne_0 = 0;

    double PFB = (x/b_FB)*(x/b_FB) + (y/b_FB)*(y/b_FB) + ((fabs(z)-z_FB)/a_FB)*((fabs(z)-z_FB)/a_FB);

    if (rLB(x, y, z) < R_LB) { // inside local bubble
        Ne_0 = J_LB * Ne_thick_disk + max(Ne_thin_disk, Ne_spiral_arms);
        if (Ne_Local_Bubbles > Ne_0) {
            return Ne_Local_Bubbles;
        } else { 
            return Ne_0;
        }
    } else if (PFB < 1.) { // inside Fermi bubble  
        Ne_0 = J_FB * Ne_thick_disk + max(Ne_thin_disk, Ne_spiral_arms);
    } else { 
        Ne_0 = Ne_thick_disk + max(Ne_thin_disk, Ne_spiral_arms);
    }

    double Ne_Gum_nebula = getNeGumNebula(pos); 
    double Ne_Loop_I = getNeLoopI(pos);

    if ((Ne_Local_Bubbles > Ne_0) 
            and (Ne_Local_Bubbles > Ne_Gum_nebula)) {
        return Ne_Local_Bubbles;
    } else if (Ne_Gum_nebula > Ne_0) {
        return Ne_Gum_nebula;
    } else if (Ne_Loop_I > Ne_0) {
        return Ne_Loop_I;
    } else {
        return Ne_0 + getNeGalacticCenter(pos);
    }
    
}

// Galactic arms ==============================================
void YMW17::generateArms() {
    // Build the arms by computing the radius of each arm with a phi step of 0.5 deg
    // the table is needed to compute the shortest distance to the arm Smin (see SjX method below)
    for (int j_arm=0; j_arm < Narms; j_arm++) {
        for (int n_phi=0; n_phi < N0_phi; n_phi++) { // step of 0.5 deg
            double phi =  n_phi * phi_step_arm + phi0_arm[j_arm]; 
            double r_j = R0_arm[j_arm] * exp((phi-phi0_arm[j_arm]) * tan(psi0_arm[j_arm]));
            radius_arm[j_arm][n_phi] = r_j;

            double nj = n0_arm[j_arm];
            if (j_arm == 2) { // Carina arm
                if (phi < phiCN) {
                    double phi_CN = (phi-phiCN) / delta_phiCN;
                    nj = nj * (1 + nCN * exp(-phi_CN*phi_CN));
                } else {
                    nj = nj * (1 + nCN);
                }
                double phi_SG = (phi-phiSG) / delta_phiSG;
                nj = nj * (1 - nSG * exp(-phi_SG*phi_SG));
            }
            n_arm[j_arm][n_phi] = nj;
        }
    }
}

double YMW17::getNeSpiralArms(const Vector3d& pos) const {
    double x =  pos.getY(); // rotate axis
    double y = -pos.getX();
    double z =  pos.getZ();
    double r =  sqrt(x*x + y*y);
    double l =  atan2(y,x); // rotate axis
    double Ne_arms = 0;

    for (int j_arm=0; j_arm < Narms; j_arm++) {
        // distance to the point (x,y) and minimal distance
        double s_min = std::numeric_limits<double>::max(); 
        double n_arm_min = 0; 
        for (int n_phi=0; n_phi < N0_phi; n_phi++) { // step of 0.5 deg
            double phi =  n_phi * phi_step_arm + phi0_arm[j_arm];
            double r_j = radius_arm[j_arm][n_phi];
            
            double x_j = -r_j * cos(phi) - x; // directly distance in x
            double y_j = -r_j * sin(phi) - y; // directly distance in y
            double s_j = sqrt(x_j*x_j + y_j*y_j); // distance point to point

            if (s_j < s_min) {
                s_min = s_j;
                n_arm_min = n_arm[j_arm][n_phi];
            }
        }

        Ne_arms = Ne_arms + n_arm_min * sech2(s_min / w_arm[j_arm]);
    }

    double H = Hfunction(r);
    return Ne_arms * sech2((z-zw(r,l))/(Ka*H)) * gd(r) * sech2((r-B2)/Aa);
}

double YMW17::getNeThickDisk(const Vector3d& pos) const {
    double x =  pos.getY(); // rotate axis
    double y = -pos.getX();
    double r = sqrt(x*x + y*y);
    double phi = atan2(y,x); // rotate axis
    double z =  pos.getZ() - zw(r,phi);
    return n1 * gd(r) * sech2(z/H1);
}

double YMW17::getNeThinDisk(const Vector3d& pos) const {
    double x =  pos.getY(); // rotate axis
    double y = -pos.getX();
    double r = sqrt(x*x + y*y);
    double phi = atan2(y,x); // rotate axis
    double z =  pos.getZ() - zw(r,phi);
    double H = Hfunction(r);
    return n2 * gd(r) * sech2((r-B2)/A2) * sech2(z/(K2*H));
}

double YMW17::Hfunction(const double R) const {
    return 32*pc + 1.6e-3 * R + 4e-7 * R*R / pc;
}

double YMW17::zw(const double R, const double phi) const {
    // Warp starts outside R_w radius
    return R < R_w ? 0 : gamma_w * (R - R_w) * cos(phi - phi_w);
}

// Galactic center ===========================================================
double YMW17::getNeGalacticCenter(const Vector3d& pos) const {
    double x =  pos.getY() - x_GC; // reverse axis
    double y = -pos.getX() - y_GC;
    double z = pos.getZ() - z_GC;
    double R_GC2 = x*x + y*y;
    return n_GC * exp(-R_GC2/(A_GC*A_GC)) * sech2((z-z_GC)/H_GC);
}

// Local Bubbles ============================================================
double YMW17::getNeLocalBubbles(const Vector3d& pos) const {
    return 0;
    //double x =  pos.getY(); // reverse axis
    //double y = -pos.getX();
    //double z =  pos.getZ();
    //double rlb = rLB(x, y, z);
    //double l = atan2(y,x); // rotate axis

    //double nLB1, nLB2;
    //// compute n_LB1
    //nLB1 = n_LB1 * sech2((l-l_LB1)/delta_l_LB1) * sech2((rlb-R_LB)/W_LB1);
    //nLB1 = nLB1 * sech2(z/H_LB1);
    //// compute n_LB2
    //nLB2 = n_LB2 * sech2((l-l_LB2)/delta_l_LB2) * sech2((rlb-R_LB)/W_LB2);
    //nLB2 = nLB2 * sech2(z/H_LB2);
    //
    //return nLB1 + nLB2; 
}

double YMW17::rLB(const double x, const double y, const double z) const {
    double R = 0.94*(y - Rsun - 40*pc) - 0.34 * z;
    return sqrt(R*R + x*x);
}

// GUM nebula ===============================================================
double YMW17::getNeGumNebula(const Vector3d& pos) const {
    return 0;
    //double x =  pos.getY() - x_GN; // reverse axis
    //double y = -pos.getX() - y_GN;
    //double z =  (pos.getZ() - z_GN) / K_GN;
    //double sn2 = (x*x + y*y + z*z) - A_GN_2;
    //return n_GN * exp(-sn2/W_GN_2);
}

// Loop I ===================================================================
double YMW17::getNeLoopI(const Vector3d& pos) const {
    return 0;
    //double x =  pos.getY() - x_LI; // reverse axis
    //double y = -pos.getX() - y_LI;
    //double z =  pos.getZ() - z_LI;
    //double rLI = sqrt(x*x + y*y + z*z);
    ////double theta = acos((x*cos_l_LI + z*sin_l_LI)/rLI) / delta_l_LI;
    //double theta = l_LI / delta_l_LI;
    //rLI = (rLI - R_LI) / W_LI;

    //return n_LI * exp(-rLI*rLI + theta*theta);
}

// useful functions ========================================================
double YMW17::sech2(const double x) const {
    double cosh_x = (exp(x) + exp(-x)) /2.; 
    return 1. / (cosh_x * cosh_x);
}

double YMW17::gd(const double R) const {
    return R < Bd ? 1 : sech2((R-Bd)/Ad);
}

} // CRPROPA NAMESPACE
