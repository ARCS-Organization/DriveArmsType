//
// Created by User on 7/10/2025.
//

#include "ArmDynamics.h"
#include <Eigen/Dense>

ArmDynamics::ArmDynamics(std::vector<double> lengthsMeters, std::vector<double> comLengthsMeters, std::vector<double> massesKilograms) {
    l1M = lengthsMeters[0];
    l2M = lengthsMeters[1];
    l3M = lengthsMeters[2];
    l4M = lengthsMeters[3];

    comL1M = comLengthsMeters[0];
    comL2M = comLengthsMeters[1];
    comL3M = comLengthsMeters[2];
    comL4M = comLengthsMeters[3];

    m1Kg = massesKilograms[0];
    m2Kg = massesKilograms[1];
    m3Kg = massesKilograms[2];
    m4Kg = massesKilograms[3];
}

Eigen::Matrix4d ArmDynamics::M(Eigen::Vector4d positionsR) {
    const double theta1 = positionsR(0);
    const double theta2 = positionsR(1);
    const double theta3 = positionsR(2);
    const double phi = positionsR(3);

    const double c12 = cos(theta1 - theta2);
    const double c13 = cos(theta1 - theta3);
    const double cphi = cos(phi - theta1 + theta3);

    Eigen::Matrix4d M(4, 4);
    M.setZero();
    // --- Mass Matrix M ---
    M << l1M*l1M*(m2Kg + m3Kg + m4Kg) + l1M*l2M*(m2Kg + m3Kg + m4Kg)*c12 + l1M*l3M*m3Kg*c13 + l1M*l4M*m4Kg*cphi + m1Kg*comL1M*comL1M,
         l1M*l2M*(m2Kg + m3Kg + m4Kg)*c12 + l1M*l1M*(m2Kg + m3Kg + m4Kg),
         l1M*l3M*m3Kg*c13 + l1M*l4M*m4Kg*cphi,
         l1M*l4M*m4Kg*cphi,

         l1M*l2M*(m2Kg + m3Kg + m4Kg)*c12 + l1M*l1M*(m2Kg + m3Kg + m4Kg),
         l1M*l2M*(m2Kg + m3Kg + m4Kg)*c12 + l2M*l2M*(m2Kg + m3Kg + m4Kg) + m2Kg*comL2M*comL2M,
         0,
         0,

         l1M*l3M*m3Kg*c13 + l1M*l4M*m4Kg*cphi,
         0,
         l1M*l3M*m3Kg*c13 + l1M*l4M*m4Kg*cphi + l3M*l3M*m3Kg + l4M*l4M*m4Kg + m3Kg*comL3M*comL3M + m4Kg*comL4M*comL4M,
         l1M*l4M*m4Kg*cphi,

         l1M*l4M*m4Kg*cphi,
         0,
         l1M*l4M*m4Kg*cphi,
         l1M*l4M*m4Kg*cphi + l4M*l4M*m4Kg + m4Kg*comL4M*comL4M;

    return M;
}

Eigen::Matrix4d ArmDynamics::C(Eigen::Vector4d positionsR, Eigen::Vector4d velocitiesRPS) {
    const double theta1 = positionsR(0);
    const double theta2 = positionsR(1);
    const double theta3 = positionsR(2);
    const double phi = positionsR(3);

    const double dtheta1 = velocitiesRPS(0);
    const double dtheta2 = velocitiesRPS(1);
    const double dtheta3 = velocitiesRPS(2);
    const double dphi = velocitiesRPS(3);

    const double s12 = sin(theta1 - theta2);
    const double s13 = sin(theta1 - theta3);
    const double sphi = sin(phi - theta1 + theta3);

    Eigen::Matrix4d C(4, 4);
    C.setZero(); // Initialize all elements to 0

    C(0, 0) =  0.5*l1M*l2M*(m2Kg + m3Kg + m4Kg)*s12*dtheta2
            - 0.5*l1M*l4M*m4Kg*sphi*dphi
            + (0.5*l1M*l3M*m3Kg*s13 - 0.5*l1M*l4M*m4Kg*sphi)*dtheta3
            + (-0.5*l1M*l2M*(m2Kg + m3Kg + m4Kg)*s12 - 0.5*l1M*l3M*m3Kg*s13 + 0.5*l1M*l4M*m4Kg*sphi)*dtheta1;

    C(0, 1) = 0.5*l1M*l2M*(m2Kg + m3Kg + m4Kg)*s12*dtheta1 + 1.5*l1M*l2M*(m2Kg + m3Kg + m4Kg)*s12*dtheta2;
    C(0, 2) = -1.5*l1M*l4M*m4Kg*sphi*dphi + (0.5*l1M*l3M*m3Kg*s13 - 0.5*l1M*l4M*m4Kg*sphi)*dtheta1
              + (1.5*l1M*l3M*m3Kg*s13 - 1.5*l1M*l4M*m4Kg*sphi)*dtheta3;
    C(0, 3) = -1.5*l1M*l4M*m4Kg*sphi*dphi - 0.5*l1M*l4M*m4Kg*sphi*dtheta1 - 1.5*l1M*l4M*m4Kg*sphi*dtheta3;

    C(1, 0) = -1.5*l1M*l2M*(m2Kg + m3Kg + m4Kg)*s12*dtheta1 - 0.5*l1M*l2M*(m2Kg + m3Kg + m4Kg)*s12*dtheta2;
    C(1, 1) = -0.5*l1M*l2M*(m2Kg + m3Kg + m4Kg)*s12*dtheta1 + 0.5*l1M*l2M*(m2Kg + m3Kg + m4Kg)*s12*dtheta2;

    C(2, 0) = 0.5*l1M*l4M*m4Kg*sphi*dphi + (-1.5*l1M*l3M*m3Kg*s13 + 1.5*l1M*l4M*m4Kg*sphi)*dtheta1
              + (-0.5*l1M*l3M*m3Kg*s13 + 0.5*l1M*l4M*m4Kg*sphi)*dtheta3;
    C(2, 2) = -0.5*l1M*l4M*m4Kg*sphi*dphi + (-0.5*l1M*l3M*m3Kg*s13 + 0.5*l1M*l4M*m4Kg*sphi)*dtheta1
              + (0.5*l1M*l3M*m3Kg*s13 - 0.5*l1M*l4M*m4Kg*sphi)*dtheta3;
    C(2, 3) = -0.5*l1M*l4M*m4Kg*sphi*dphi + 0.5*l1M*l4M*m4Kg*sphi*dtheta1 - 0.5*l1M*l4M*m4Kg*sphi*dtheta3;

    C(3, 0) = 0.5*l1M*l4M*m4Kg*sphi*dphi + 1.5*l1M*l4M*m4Kg*sphi*dtheta1 + 0.5*l1M*l4M*m4Kg*sphi*dtheta3;
    C(3, 2) = -0.5*l1M*l4M*m4Kg*sphi*dphi + 0.5*l1M*l4M*m4Kg*sphi*dtheta1 - 0.5*l1M*l4M*m4Kg*sphi*dtheta3;
    C(3, 3) = -0.5*l1M*l4M*m4Kg*sphi*dphi + 0.5*l1M*l4M*m4Kg*sphi*dtheta1 - 0.5*l1M*l4M*m4Kg*sphi*dtheta3;

    return C;
}

Eigen::Vector4d ArmDynamics::G(Eigen::Vector4d positionsR) {
    const double theta1 = positionsR(0);
    const double theta2 = positionsR(1);
    const double theta3 = positionsR(2);
    const double phi = positionsR(3);

    Eigen::Vector4d G(4);
    G.setZero();
    G(0) = g*l1M*(m2Kg + m3Kg + m4Kg)*cos(theta1) + g*m1Kg*comL1M*cos(theta1);
    G(1) = g*l2M*(m2Kg + m3Kg + m4Kg)*cos(theta2) + g*m2Kg*comL2M*cos(theta2);
    G(2) = g*l3M*m3Kg*cos(theta3) + g*l4M*m4Kg*cos(phi + theta3)
         + g*m3Kg*comL3M*cos(theta3) + g*m4Kg*comL4M*cos(phi + theta3);
    G(3) = g*l4M*m4Kg*cos(phi + theta3) + g*m4Kg*comL4M*cos(phi + theta3);

    return G;
}

Eigen::Vector4d ArmDynamics::getAccelerationsRPSS(
        Eigen::Vector4d torqueNm,
        Eigen::Vector4d positionsR,
        Eigen::Vector4d velocitiesRPS) {
    return
        ((M(positionsR).inverse()))
            * (torqueNm - (C(positionsR, velocitiesRPS) * velocitiesRPS) - G(positionsR));
};

Eigen::Vector4d ArmDynamics::getTorquesNm(
        Eigen::Vector4d accelerationRPSS,
        Eigen::Vector4d positionsR,
        Eigen::Vector4d velocitiesRPS) {
    return
        (M(positionsR) * accelerationRPSS)
        + (C(positionsR, velocitiesRPS) * velocitiesRPS)
        + G(positionsR);
};