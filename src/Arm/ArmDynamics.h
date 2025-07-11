//
// Created by User on 7/10/2025.
//

#ifndef ARMDYNAMICS_H
#define ARMDYNAMICS_H


#include <Eigen/Dense>
#include <cmath>

class ArmDynamics {
private:
    // Generalized coordinates
    double theta1, theta2, theta3, phi;
    double dtheta1, dtheta2, dtheta3, dphi;

    // Link lengths and center of mass distances
    double l1, l2, l3, l4;
    double lambda1, lambda2, lambda3, lambda4;

    // Masses
    double m1, m2, m3, m4;

    // Gravity
    double g;

    // Matrices
    Eigen::Matrix4d M;
    Eigen::Matrix4d C;
    Eigen::Vector4d G;

public:
    ArmDynamics(); // Constructor to initialize variables

    void computeDynamics(); // Populates M, C, G
    Eigen::Matrix4d getMassMatrix() const { return M; }
    Eigen::Matrix4d getCoriolisMatrix() const { return C; }
    Eigen::Vector4d getGravityVector() const { return G; }
};



#endif //ARMDYNAMICS_H
