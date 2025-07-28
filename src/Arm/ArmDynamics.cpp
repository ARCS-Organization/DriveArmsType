//
// Created by User on 7/10/2025.
//

#include "ArmDynamics.h"

#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

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
            * (
            torqueNm
            -
            C(positionsR, velocitiesRPS) * velocitiesRPS
            -
            G(positionsR));
}

Eigen::Vector4d ArmDynamics::getNoTorquesAccelerationRPSS(Eigen::Vector4d positionsR, Eigen::Vector4d velocitiesRPS) {
    return
        ((M(positionsR).inverse())) * (
            -
            C(positionsR, velocitiesRPS) * velocitiesRPS
            -
            G(positionsR));
}

Eigen::Matrix4d ArmDynamics::noTorquesJacobianPosition(Eigen::Vector4d positionsR, Eigen::Vector4d velocitiesRPS) {
    double epsilon = 0.001;

    Eigen::Vector4d vec1 = Eigen::Vector4d();
    Eigen::Vector4d newPos1 = positionsR;

    Eigen::Vector4d vec2 = Eigen::Vector4d();
    Eigen::Vector4d newPos2 = positionsR;

    Eigen::Vector4d delta = Eigen::Vector4d();

    Eigen::Matrix4d dE = Eigen::Matrix4d();

    for(int i = 0; i < 4; i++) {
        newPos1 = positionsR;
        newPos1(i) = positionsR(i) - epsilon;
        vec1 = getNoTorquesAccelerationRPSS(newPos1, velocitiesRPS);

        newPos2 = positionsR;
        newPos2(i) = positionsR(i) + epsilon;
        vec2 = getNoTorquesAccelerationRPSS(newPos2, velocitiesRPS);

        delta = (vec2 - vec1) / epsilon;

        for(int j = 0; j < 4; j++) dE(j, i) = delta(i);
    }

    return dE;
}

Eigen::Matrix4d ArmDynamics::noTorquesJacobianVelocity(Eigen::Vector4d positionsR, Eigen::Vector4d velocitiesRPS) {
    double epsilon = 0.001;

    Eigen::Vector4d vec1 = Eigen::Vector4d();
    Eigen::Vector4d newVel1 = velocitiesRPS;

    Eigen::Vector4d vec2 = Eigen::Vector4d();
    Eigen::Vector4d newVel2 = velocitiesRPS;

    Eigen::Vector4d delta = Eigen::Vector4d();

    Eigen::Matrix4d dE = Eigen::Matrix4d();

    for(int i = 0; i < 4; i++) {
        newVel1 = velocitiesRPS;
        newVel1(i) = velocitiesRPS(i) - epsilon;
        vec1 = getNoTorquesAccelerationRPSS(positionsR, newVel1);

        newVel2 = velocitiesRPS;
        newVel2(i) = velocitiesRPS(i) + epsilon;
        vec2 = getNoTorquesAccelerationRPSS(positionsR, newVel2);

        delta = (vec2 - vec1) / epsilon;

        for(int j = 0; j < 4; j++) dE(j, i) = delta(i);
    }

    return dE;
}

Eigen::Vector4d ArmDynamics::getTorquesNm(
        Eigen::Vector4d accelerationRPSS,
        Eigen::Vector4d positionsR,
        Eigen::Vector4d velocitiesRPS) {
    return
        (M(positionsR) * accelerationRPSS)
        + (C(positionsR, velocitiesRPS) * velocitiesRPS)
        + G(positionsR);
}

Eigen::Vector<double, 8> ArmDynamics::x(Eigen::Vector4d positionsR, Eigen::Vector4d velocitiesRPSS) {
    Eigen::Vector<double, 8> x = Eigen::VectorXd();
    for(int i = 0; i < 4; i++) x(i) = positionsR(i);
    for (int i = 0; i < 4; i++) x(i+4) = velocitiesRPSS(i);
    return x;
}

Eigen::Matrix<double, 8, 8> ArmDynamics::A(Eigen::Vector4d positionsR, Eigen::Vector4d velocitiesRPS) {
    Eigen::Vector4d oldPositionR = positionsR;
    Eigen::Vector4d oldVelocityRPSS = velocitiesRPS;

    Eigen::Matrix<double, 8, 8> A = Eigen::Matrix<double, 8, 8>();
    A.setZero();

    for(int i = 0; i < 4; i++) A(i, 4+i) = 1;

    double epsilon = 0.2;

    Eigen::Vector4d vec1;
    Eigen::Vector4d newPos1;
    Eigen::Vector4d newVel1;

    Eigen::Vector4d vec2;
    Eigen::Vector4d newPos2;
    Eigen::Vector4d newVel2;

    Eigen::Vector4d delta;

    for(int i = 0; i < 4; i++) {
        newPos1 = positionsR;
        newPos1(i) = positionsR(i) - epsilon;
        vec1 = getNoTorquesAccelerationRPSS(newPos1, velocitiesRPS);

        newPos2 = positionsR;
        newPos2(i) = positionsR(i) + epsilon;
        vec2 = getNoTorquesAccelerationRPSS(newPos2, velocitiesRPS);

        delta = (vec2 - vec1) / epsilon;

        for(int j = 0; j < 4; j++)  A(4+j, i) = delta(j);

        newVel1 = velocitiesRPS;
        newVel1(i) = velocitiesRPS(i) - epsilon;
        vec1 = getNoTorquesAccelerationRPSS(positionsR, newVel1);

        newVel2 = velocitiesRPS;
        newVel2(i) = oldVelocityRPSS(i) + epsilon;
        vec2 = getNoTorquesAccelerationRPSS(positionsR, newVel2);

        delta = (vec2 - vec1) / epsilon;

        for(int j = 0; j < 4; j++) A(4+j, 4+i) = delta(j);
    }

    return A;
}

Eigen::Matrix<double, 8, 8> ArmDynamics::B(Eigen::Vector4d positionsR) {
    Eigen::Matrix<double, 8, 8> B = Eigen::Matrix<double, 8, 8>();
    B.setZero();
    Eigen::Matrix4d mInv = M(positionsR).inverse();
    for(int i = 0; i < 4; i++) for(int j = 0; j < 4; j++) B(i+4, j) = mInv(i, j);
    return B;
}

Eigen::Matrix<double, 4, 8> ArmDynamics::C() {
    Eigen::Matrix<double, 4, 8> C = Eigen::Matrix<double, 4, 8>();
    C.setZero();
    for(int i = 0; i < 4; i++) C(i,i) = 1;
    return C;
}

Eigen::Vector4d ArmDynamics::D() {
    Eigen::Vector4d D = Eigen::Vector4d();
    D.setZero();
    return D;
}

Eigen::Matrix<double, 8, 8> ArmDynamics::discretizedA(Eigen::Vector4d positionsR, Eigen::Vector4d velocityRPSS, double dt) {
    return (ArmDynamics::A(positionsR, velocityRPSS) * dt).exp();
}

Eigen::Matrix<double, 8, 8> ArmDynamics::discretizedB(Eigen::Vector4d positionsR, Eigen::Matrix<double, 8, 8> Ac, Eigen::Matrix<double, 8, 8> Ad, double dt) {
    Eigen::Matrix<double, 8, 8> I = Eigen::Matrix<double, 8, 8>();
    I.setIdentity();

    std::cout << Ac << std::endl;
    std::cout << Ac.determinant() << std::endl;
    std::cout << Ac.inverse().determinant() << std::endl;
    std::cout << Ad - I << std::endl;
    std::cout << B(positionsR) << std::endl;


    return Ac.inverse() * (Ad - I) * (B(positionsR));
}