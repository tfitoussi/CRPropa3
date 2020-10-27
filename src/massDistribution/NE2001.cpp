#include "crpropa/massDistribution/NE2001.h"
#include <algorithm>
#include <fstream>

namespace crpropa {
using namespace std;

NE2001::NE2001(){
    // Common parameters
    Rsun = 8.5 * kpc;
    
    // Parameters Thick Disk
    //H1 = 0.97 * kpc; 
    //n1 = 0.033*kpc / H1; // cm^-3
    H1 = 1.31 * kpc;  // Correction from Schnitzeler 2012
    n1 = 0.016*kpc / H1; // cm^-3
    A1 = 17.5 * kpc;

    // Parameters Thin Disk
    n2 = 0.08; // cm^-3
    H2 = 0.15 * kpc;
    A2 = 3.8 * kpc;
    
    // Parameters Spiral Arms
    Narms = 5;
    Aa = 10.5 * kpc;
    na = 0.028;
    n_arm[0] = 0.50 * na; // cm^-3 
    n_arm[1] = 1.20 * na; // cm^-3 
    n_arm[2] = 1.30 * na; // cm^-3 
    n_arm[3] = 1.00 * na; // cm^-3 
    n_arm[4] = 0.25 * na; // cm^-3
    wa = 0.65 * kpc;
    w_arm[0] = 1.0 * wa;
    w_arm[1] = 1.5 * wa;
    w_arm[2] = 1.0 * wa;
    w_arm[3] = 0.8 * wa;
    w_arm[4] = 1.0 * wa;
    Ha = 0.23 * kpc;
    h_arm[0] = 1.0 * Ha;
    h_arm[1] = 0.8 * Ha;
    h_arm[2] = 1.3 * Ha;
    h_arm[3] = 1.5 * Ha;
    h_arm[4] = 1.0 * Ha;
    alpha_arm[0] = 4.25;
    alpha_arm[1] = 4.25;
    alpha_arm[2] = 4.89;
    alpha_arm[3] = 4.89;
    alpha_arm[4] = 4.57;
    R_min_arm[0] = 3.48 * kpc;
    R_min_arm[1] = 3.48 * kpc;
    R_min_arm[2] = 4.90 * kpc;
    R_min_arm[3] = 3.76 * kpc;
    R_min_arm[4] = 8.10 * kpc;
    theta_min_arm[0] = 0.000; // rad
    theta_min_arm[1] = 3.141; // rad
    theta_min_arm[2] = 2.525; // rad
    theta_min_arm[3] = 4.240; // rad
    theta_min_arm[4] = 5.847; // rad
    extent_arm[0] = 6.00; // rad
    extent_arm[1] = 6.00; // rad
    extent_arm[2] = 6.00; // rad
    extent_arm[3] = 6.00; // rad
    extent_arm[4] = 0.55; // rad
    // to deform arm 2 and 4
    theta_arm2_lim[0] = 180. * deg; // rad 
    theta_arm2_lim[1] = 315. * deg; // rad  
    theta_arm2_lim[2] = 370. * deg; // rad  
    theta_arm2_lim[3] = 410. * deg; // rad  
    theta_arm2[0] = 260. * deg; // rad  
    theta_arm2[1] = 345. * deg; // rad  
    theta_arm2[2] = 390. * deg; // rad  
    theta_arm4_lim[0] = 290. * deg; // rad  
    theta_arm4_lim[1] = 395. * deg; // rad  
    theta_arm4 = 350. * deg; // rad   
    theta2a = 340. * deg; // rad   
    theta2b = 370. * deg; // rad    
    theta3a = 290. * deg; // rad   
    theta3b = 363. * deg; // rad    
    theta_step_arm = 0.5*deg;
    generateArms();

    // Parameters Galactic Center (GC)
    n_gc0 = 10.0; // cm^-3 
    R_gc = 0.145 * kpc;
    H_gc = 0.26 * kpc;
    x_gc = -0.01 * kpc; 
    y_gc = 0. * kpc; 
    z_gc = -0.02 * kpc;
    
    // Parameters Voids
    loadNeVoid(); // load voids parameters 
    
    // Parameters Clumps
    loadNeClumps(); // load voids parameters 
    
    // Parameters Local Interstellar Medium (LISM)
    double aldr = 1.50 * kpc; 
    double bldr = 0.75 * kpc; 
    double cldr = 0.50 * kpc;
    double thetaldr = -24.2 * deg; 
    double cos_thetaldr = cos(thetaldr);
    double sin_thetaldr = sin(thetaldr);
    apdr = cos_thetaldr*cos_thetaldr/aldr/aldr + sin_thetaldr*sin_thetaldr/bldr/bldr;
    bpdr = sin_thetaldr*sin_thetaldr/aldr/aldr + cos_thetaldr*cos_thetaldr/bldr/bldr;
    cpdr = 1. / cldr / cldr;
    dpdr = 2 * sin_thetaldr * cos_thetaldr * (1./aldr/aldr - 1./bldr/bldr);
    xldr = 1.36 * kpc; 
    yldr = 8.06 * kpc; 
    zldr = 0.00 * kpc;
    neldr0 = 0.0120; 

    double alsb = 1.05 * kpc; 
    double blsb = 0.425 * kpc; 
    double clsb = 0.325 * kpc;
    double thetalsb = 139. * deg; 
    double cos_thetalsb = cos(thetalsb);
    double sin_thetalsb = sin(thetalsb);
    apsb = cos_thetalsb*cos_thetalsb/alsb/alsb + sin_thetalsb*sin_thetalsb/blsb/blsb;
    bpsb = sin_thetalsb*sin_thetalsb/alsb/alsb + cos_thetalsb*cos_thetalsb/blsb/blsb;
    cpsb = 1. / clsb / clsb;
    dpsb = 2 * sin_thetalsb * cos_thetalsb * (1./alsb/alsb - 1./blsb/blsb);
    xlsb = -0.75 * kpc; 
    ylsb = 9.0 * kpc; 
    zlsb = -0.05 * kpc;
    nelsb0 = 0.016; 

    alhb = 0.085 * kpc; 
    blhb = 0.1 * kpc; 
    clhb = 0.33 * kpc; 
    double thetalhb = 15 * deg;
    yzslope = tan(thetalhb);
    xlhb = 0.01 * kpc; 
    ylhb = 8.45 * kpc; 
    zlhb = 0.17 * kpc;
    nelhb0 = 0.005; 

    xlpI = -0.045 * kpc; 
    ylpI = 8.4 * kpc; 
    zlpI = 0.07 * kpc;
    rlpI = 0.12; 
    double drlpI = 0.06;
    alpI = rlpI + drlpI;
    nelpI = 0.0125; 
    dnelpI = 0.0125; 
}

double NE2001::getNe(const Vector3d& pos) const  {
    double Ne;
    double Ne_voids = getNeVoids(pos);

    if (Ne_voids > 0) { // inside a void (w_void = 0)
        Ne = Ne_voids;
    } else {
        double Ne_lism = getNeLISM(pos); 
        if (Ne_lism > 0) { // inside Local interstellar medium
            Ne = Ne_lism;
        } else {
            Ne = getNeDisk(pos) + getNeGalacticCenter(pos);
        }
    }
    Ne = Ne + getNeClumps(pos); 
    return Ne;
}

// Voids =======================================================
void NE2001::loadNeVoid() {
    double l_void;
    double b_void;
    double d_void;
    double F_void;
    double th1; 
    double th2; 
    int flag, edge; 
    string name;

    // load parameters 
    ifstream datafile;
    string line;
    datafile.open(getDataPath("NE2001/nevoidN.dat"));
    // remove first line
    getline(datafile, line);

    N_voids = 0;
    while (datafile >> flag >> l_void >> b_void >> d_void >> ne_void[N_voids] >> F_void >> aa_void[N_voids] >> bb_void[N_voids] >> cc_void[N_voids] >> th1 >> th2 >> edge >> name) {
        // rotation matrix in the 'q = ' statement getNeVoids
        // rotation matrix corresponds to \Lambda_z\Lambda_y
        // where \Lambda_y = rotation around y axis
        //       \Lambda_z = rotation around z axis
        // defined as
        // \Lambda_y =  c1  0  s1 
        //               0  1   0
        //             -s1  0  c1
        //
        // \Lambda_z =  c2 s2   0 
        //             -s2 c2   0
        //               0  0   1
        // =>
        // \Lambda_z\Lambda_y =  c1*c2   s2   s1*c2
        //                      -s2*c1   c2  -s1*s2
        //                         -s1    0      c1  
        //
        // so the rotation is around the y axis first, then the z axis

        aa_void[N_voids] = aa_void[N_voids] * kpc;
        bb_void[N_voids] = bb_void[N_voids] * kpc;
        cc_void[N_voids] = cc_void[N_voids] * kpc;

        d_void = d_void * kpc;
        l_void = l_void * deg; // convert in rad
        b_void = b_void * deg; // convert in rad
        x_void[N_voids] = d_void * cos(b_void) * sin(l_void);
        y_void[N_voids] = Rsun - d_void * cos(b_void) * cos(l_void);
        z_void[N_voids] = d_void * sin(b_void);

        th1 = th1 * deg; // convert in rad
        th2 = th2 * deg; // convert in rad
        c2_void[N_voids] = cos(th2);
        s2_void[N_voids] = sin(th2);
        c1_void[N_voids] = cos(th1);
        s1_void[N_voids] = sin(th1);
        c1c2_void[N_voids] = c1_void[N_voids] * c2_void[N_voids];
        s1c2_void[N_voids] = s1_void[N_voids] * c2_void[N_voids];
        c1s2_void[N_voids] = c1_void[N_voids] * s2_void[N_voids];
        s1s2_void[N_voids] = s1_void[N_voids] * s2_void[N_voids];

        N_voids++;
    }
    datafile.close();
}

double NE2001::getNeVoids(const Vector3d& pos) const {
    double x = pos.getY(); // reverse axis
    double y = -pos.getX();
    double z = pos.getZ();
    double Nev = 0.;
    // rotation matrix define in loadNeVoids
    for (int iv=0; iv<N_voids; iv++) {
        double dx = x - x_void[iv];
        double dy = y - y_void[iv];
        double dz = z - z_void[iv];
        double term1, term2, term3;
        term1 = (c1c2_void[iv] *dx + s2_void[iv]*dy + s1c2_void[iv]*dz) / aa_void[iv];
        term2 = (-c1s2_void[iv]*dx + c2_void[iv]*dy - s1s2_void[iv]*dz) / bb_void[iv];
        term3 = (-s1_void[iv]*dx + c1_void[iv]*dz) / cc_void[iv];
        double q = term1*term1 + term2*term2 + term3*term3;

        if (q<=1.) {
            Nev = Nev + ne_void[iv];
        }
    }
    return Nev *cm*cm*cm;
}

// Galactic arms ==============================================
void NE2001::generateArms() {
    // Build the arms by computing the radius of each arm with a theta step of 0.5 deg
    // the table is needed to compute the shortest distance to the arm Smin (see SjX method below)
    for (int j_arm=0; j_arm < Narms; j_arm++) {
        double theta = theta_min_arm[j_arm];
        n_theta_arm[j_arm] = 0;
        while (theta <= theta_min_arm[j_arm] + extent_arm[j_arm]) {
            double radius;
            radius = R_min_arm[j_arm] * exp((theta - theta_min_arm[j_arm])/alpha_arm[j_arm]);
            // arm 2 and 4 are deformed  locally 
            // (arm numbered from 1 to 5 but j_arm from 0 to 4)
            if (j_arm+1 == 2) {
                if ((theta > theta_arm2_lim[0]) and (theta <= theta_arm2_lim[1])) {
                    radius = radius * (1. + 0.16 * cos((theta - theta_arm2[0])*180./135.));
                } else if ((theta > theta_arm2_lim[1]) and (theta <= theta_arm2_lim[2])) {
                    radius = radius * (1. - 0.07 * cos((theta - theta_arm2[1])*180./55.));
                } else if ((theta > theta_arm2_lim[2]) and (theta <= theta_arm2_lim[3])) {
                    radius = radius * (1. + 0.04 * cos((theta - theta_arm2[2])*180./40.));
                }
            } else if (j_arm+1 == 4) {
                if ((theta > theta_arm4_lim[0]) and (theta <= theta_arm4_lim[1])) {
                    radius = radius * (1. - 0.11 * cos((theta - theta_arm4)*180./105.));
                }
            }
            // better performance with static array than std::vector
            radius_arm[j_arm][n_theta_arm[j_arm]] = radius;
            theta += theta_step_arm;
            n_theta_arm[j_arm]++;
        }
    }
}

double NE2001::SjX(const double x, const double y, const int j_arm) const {
    double r_j, x_j, y_j; // coordinates along the arm
    double s_j, s_min = std::numeric_limits<double>::max(); // distance to the point (x,y) and minimal distance

    for (int n_theta=0; n_theta < n_theta_arm[j_arm]; n_theta++) { // step of 0.5 deg
        double theta = theta_min_arm[j_arm] + n_theta * theta_step_arm;
        double radius;
        r_j = radius_arm[j_arm][n_theta];
        x_j = -r_j * cos(theta) - x; // directly distance in x
        y_j = -r_j * sin(theta) - y; // directly distance in y
        s_j = sqrt(x_j*x_j + y_j*y_j); // distance point to point

        if (s_j < s_min) s_min = s_j;
    }
    return s_min;
}

double NE2001::getNeDisk(const Vector3d& pos) const {
    double Ne_thick_disk = getNeThickDisk(pos);
    double Ne_thin_disk = getNeThinDisk(pos);
    double Ne_spiral_arms = getNeSpiralArms(pos);
    return Ne_thick_disk + Ne_thin_disk + Ne_spiral_arms;
}

double NE2001::getNeThickDisk(const Vector3d& pos) const {
    double x =  pos.getY(); // reverse axis
    double y = -pos.getX();
    double r = sqrt(x*x + y*y);
    double g1 = cos((M_PI*r)/(2*A1)) /  cos((M_PI*Rsun)/(2*A1));
    double h1 = sech2(pos.getZ() / H1);
    return n1 * g1 * h1 *cm*cm*cm;
}

double NE2001::getNeThinDisk(const Vector3d& pos) const {
    double x =  pos.getY(); // reverse axis
    double y = -pos.getX();
    double r = sqrt(x*x + y*y);
    double r_arg = (r-A2) / (1.8*kpc);
    double g2 = exp(- r_arg*r_arg);
    double h2 = sech2(pos.getZ() / H2);
    return n2 * g2 * h2 *cm*cm*cm;
}

double NE2001::getNeSpiralArms(const Vector3d& pos) const {
    double x =  pos.getY(); // reverse axis
    double y = -pos.getX();
    double z = pos.getZ();
    double r = sqrt(x*x + y*y);
    double theta = atan2(y,x);
    double Ne_arms = 0;
    for (int j=0; j < Narms; j++) {
        double sjx = SjX(pos.x, pos.y, j)/w_arm[j];
        double g_aj = exp(-sjx*sjx);
        if (r > Aa) g_aj = g_aj * sech2((r-Aa)/(2.*kpc));
        if ((j+1 == 2) and (theta >= theta2a) and (theta <= theta2b)) {
            double arg = 2.*M_PI *(theta - theta2a)/(theta2b - theta2a);
            double fac = (1.1 + 0.9 * cos(arg)) / 2.;
            g_aj = g_aj * fac;
        } else if ((j+1 == 3) and (theta >= theta3a) and (theta <= theta3b)) {
            double arg = 2.*M_PI *(theta - theta3a)/(theta3b - theta3a);
            double fac = (1. + cos(arg)) / 2.;
            g_aj = g_aj * fac * fac * fac * fac; // g_aj * fac**4 
        }
        Ne_arms = Ne_arms + n_arm[j] * g_aj * sech2(z / h_arm[j]); 
    }
    return Ne_arms *cm*cm*cm;
}

double NE2001::sech2(const double x) const {
    double cosh_x = (exp(x) + exp(-x)) /2.; 
    return 1. / (cosh_x * cosh_x);
}

// Galactic center ===========================================================
double NE2001::getNeGalacticCenter(const Vector3d& pos) const {
    double x =  pos.getY() - x_gc; // reverse axis
    double y = -pos.getX() - y_gc;
    double z = pos.getZ() - z_gc;
    double delta_Rperp2 = x*x + y*y;
    double ne_gc = n_gc0 * exp(-delta_Rperp2/(R_gc*R_gc));
    ne_gc = ne_gc * exp(-z*z/(H_gc*H_gc)); 
    return ne_gc *cm*cm*cm;
}

// Clumps ===================================================================
void NE2001::loadNeClumps() {
    // LOS quantities. 
    // lc,bc = Galactic coordinates (deg)
    //   nec = clump electron density (cm^{-3})
    //    Fc = fluctuation parameter
    //    dc = clump distance from Earth (kpc)
    //    rc = clump radius (kpc)
    //  edge = 0,1  0=> Gaussian, 1=> Gaussian w/ hard edge at e^{-1} 
    //  type = LOS type (P pulsar, G other Galactic, X extragalactic
    // losname = useful name
    double l_clumps, b_clumps;
    double d_clumps;
    double F_clumps;
    int flag, edge; 
    string type, name;

    // load parameters 
    ifstream datafile;
    string line;
    datafile.open(getDataPath("NE2001/neclumpN.dat"));
    // remove first line
    getline(datafile, line);
    
    N_clumps = 0;
    while (datafile >> flag >> l_clumps >> b_clumps >> ne_clumps[N_clumps] >> F_clumps >> d_clumps >> R_clumps[N_clumps] >> edge >> type >> name) {
        R_clumps[N_clumps] = R_clumps[N_clumps] * kpc;
        d_clumps = d_clumps * kpc;
        l_clumps = l_clumps * deg; // convert in rad
        b_clumps = b_clumps * deg; // convert in rad
        x_clumps[N_clumps] = d_clumps * cos(b_clumps) * sin(l_clumps);
        y_clumps[N_clumps] = Rsun - d_clumps * cos(b_clumps) * cos(l_clumps);
        z_clumps[N_clumps] = d_clumps * sin(b_clumps);        

        N_clumps++;
    }
    datafile.close();
}

double NE2001::getNeClumps(const Vector3d& pos) const {
    double x =  pos.getY(); // reverse axis
    double y = -pos.getX();
    double z = pos.getZ();
    double Nec = 0.;
    // rotation matrix define in loadNeVoids
    for (int ic=0; ic<N_clumps; ic++) {
        double dx = (x - x_clumps[ic]) / R_clumps[ic];
        double dy = (y - y_clumps[ic]) / R_clumps[ic];
        double dz = (z - z_clumps[ic]) / R_clumps[ic];
        double q = dx*dx + dy*dy +dz*dz;
        if (q<=1.) {
            Nec = Nec + ne_clumps[ic];
        }
    }
    return Nec *cm*cm*cm;
}

// Local Interstellar Medium ================================================
double NE2001::getNeLISM(const Vector3d& pos) const {
    //
    //  Beware that w_LISM is computed inside the method ...
    //  Need to find a better way of doing it ...
    //
    // method to calculate the electron density for the 
    // Local Interstellar Medium
    //
    // JMC 26 August-11 Sep. 2000
    //     25 October 2001: modified to change weighting scheme
    //                      so that the ranking is LHB: LSB: LDR
    //                      (LHB overrides LSB and LDR; LSB overrides LDR)
    //     16 November 2001: added Loop I component with weighting scheme
    //		        LHB:LOOPI:LSB:LDR	
    //		        LHB   overides everything,
    //			LOOPI overrides LSB and LDR
    //			LSB   overrides LDR
    //			LISM  overrides general Galaxy
    //     20 November 2001: The LOOPI component is truncated below z=0 
    //
    // after discussions with Shami Chatterjee
    // the sizes, locations and densities of the LISM components
    // are based largely on work by Toscano et al. 1999
    // and Bhat et al. 1999
    double x =  pos.getY(); // reverse axis
    double y = -pos.getX();
    double z = pos.getZ();
    double q;

    // low density region in Q1
    q = (x-xldr)*(x-xldr)*apdr + (y-yldr)*(y-yldr)*bpdr + (z-zldr)*(z-zldr)*cpdr + (x-xldr)*(y-yldr)*dpdr; 
    double w_ldrq1 = q <= 1. ? 1. : 0.; 
    double  Ne_ldrq1 = q <= 1. ? neldr0 : 0.;

    // Local Super Bubble
    q = (x-xlsb)*(x-xlsb)*apsb + (y-ylsb)*(y-ylsb)*bpsb + (z-zlsb)*(z-zlsb)*cpsb + (x-xlsb)*(y-ylsb)*dpsb; 
    double w_lsb = q <= 1. ? 1. : 0.; 
    double  Ne_lsb = q <= 1. ? nelsb0 : 0.;

    // Local Hot Bubble
    double aa = alhb;
    if ((z <= 0.) & (z >= (zlhb-clhb))){
        aa = 0.001 + (alhb-0.001)*(1. - (1./(zlhb-clhb))*z);
    }
    double yaxis = ylhb + yzslope*z;  
    double qxy = (x-xlhb)*(x-xlhb)/aa/aa + (y-yaxis)*(y-yaxis)/blhb/blhb;
	double qz =  fabs(z-zlhb)/clhb;
    bool cond = (qxy <=1.) & (qz <= 1.);
    double w_lhb = cond ? 1. : 0.; 
    double  Ne_lhb = cond ? nelhb0 : 0.;

    // Loop I
    double w_loopI, Ne_loopI;
    if (z < 0){
        w_loopI = 0.;
        Ne_loopI = 0.;
    } else {
        double r = sqrt((x-xlpI)*(x-xlpI) + (y-ylpI)*(y-ylpI) + (z-zlpI)*(z-zlpI));
        if (r > alpI) {
            w_loopI = 0.;
            Ne_loopI = 0.;
        } else if (r <= rlpI) {
            w_loopI = 1.;
            Ne_loopI = nelpI;
        } else {
            w_loopI = 1.;
            Ne_loopI = dnelpI;
        }
    }

    // Total
    double Ne_LISM = 0.;
	double w_LISM = max(w_loopI, max(w_ldrq1, max(w_lsb, w_lhb)));
    if (w_LISM > 0) {
        Ne_LISM = (1-w_loopI) * (w_lsb * Ne_lsb + (1-w_lsb)*Ne_ldrq1);
        Ne_LISM = Ne_LISM + w_loopI * Ne_loopI;
        Ne_LISM = (1-w_lhb) * Ne_LISM;
        Ne_LISM = Ne_LISM + w_lhb * Ne_lhb;
    }

    return Ne_LISM *cm*cm*cm;
}

} // CRPROPA NAMESPACE
