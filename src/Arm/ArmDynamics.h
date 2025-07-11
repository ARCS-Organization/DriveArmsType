//
// Created by User on 7/10/2025.
//

#ifndef ARMDYNAMICS_H
#define ARMDYNAMICS_H


#include <Eigen/Dense>
#include <cmath>

class ArmDynamics {
private:
    // Link lengths and center of mass distances
    double l1M, l2M, l3M, l4M;
    double comL1M, comL2M, comL3M, comL4M;

    // Masses
    double m1Kg, m2Kg, m3Kg, m4Kg;

    // Gravity
    const double g = 9.81;

public:
    ArmDynamics(std::vector<double> lengthsMeters, std::vector<double> comLengthsMeters, std::vector<double> massesKilograms); // Constructor to initialize variables

    Eigen::Matrix4d M(Eigen::Vector4d positionsR);
    Eigen::Matrix4d C(Eigen::Vector4d positionsR, Eigen::Vector4d velocitiesRPS);
    Eigen::Vector4d G(Eigen::Vector4d positionsR);

    /* For simulations */
    Eigen::Vector4d getAccelerationsRPSS(
        Eigen::Vector4d torqueNm,
        Eigen::Vector4d positionsR,
        Eigen::Vector4d velocitiesRPS);

    /* For controls */
    Eigen::Vector4d getTorquesNm(
        Eigen::Vector4d accelerationRPSS,
        Eigen::Vector4d positionsR,
        Eigen::Vector4d velocitiesRPS);
};


#endif //ARMDYNAMICS_H
