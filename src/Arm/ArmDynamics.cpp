//
// Created by User on 7/10/2025.
//

#include "ArmDynamics.h"

ArmDynamics::ArmDynamics() {

}



void ArmDynamics::computeDynamics() {
    // --- Shorthand ---
    double c12 = cos(theta1 - theta2);
    double c13 = cos(theta1 - theta3);
    double cphi = cos(phi - theta1 + theta3);
    double cphiT = cos(phi + theta3);
    double s12 = sin(theta1 - theta2);
    double s13 = sin(theta1 - theta3);
    double sphi = sin(phi - theta1 + theta3);
    double sphiT = sin(phi + theta3);

    // --- Mass Matrix M ---
    M << l1*l1*(m2 + m3 + m4) + l1*l2*(m2 + m3 + m4)*c12 + l1*l3*m3*c13 + l1*l4*m4*cphi + m1*lambda1*lambda1,
         l1*l2*(m2 + m3 + m4)*c12 + l1*l1*(m2 + m3 + m4),
         l1*l3*m3*c13 + l1*l4*m4*cphi,
         l1*l4*m4*cphi,

         l1*l2*(m2 + m3 + m4)*c12 + l1*l1*(m2 + m3 + m4),
         l1*l2*(m2 + m3 + m4)*c12 + l2*l2*(m2 + m3 + m4) + m2*lambda2*lambda2,
         0, 0,

         l1*l3*m3*c13 + l1*l4*m4*cphi,
         0,
         l1*l3*m3*c13 + l1*l4*m4*cphi + l3*l3*m3 + l4*l4*m4 + m3*lambda3*lambda3 + m4*lambda4*lambda4,
         l1*l4*m4*cphi,

         l1*l4*m4*cphi,
         0,
         l1*l4*m4*cphi,
         l1*l4*m4*cphi + l4*l4*m4 + m4*lambda4*lambda4;

    // --- Coriolis Matrix C ---
    C.setZero(); // Initialize all elements to 0

    C(0, 0) =  0.5*l1*l2*(m2 + m3 + m4)*s12*dtheta2
            - 0.5*l1*l4*m4*sphi*dphi
            + (0.5*l1*l3*m3*s13 - 0.5*l1*l4*m4*sphi)*dtheta3
            + (-0.5*l1*l2*(m2 + m3 + m4)*s12 - 0.5*l1*l3*m3*s13 + 0.5*l1*l4*m4*sphi)*dtheta1;

    C(0, 1) = 0.5*l1*l2*(m2 + m3 + m4)*s12*dtheta1 + 1.5*l1*l2*(m2 + m3 + m4)*s12*dtheta2;
    C(0, 2) = -1.5*l1*l4*m4*sphi*dphi + (0.5*l1*l3*m3*s13 - 0.5*l1*l4*m4*sphi)*dtheta1
              + (1.5*l1*l3*m3*s13 - 1.5*l1*l4*m4*sphi)*dtheta3;
    C(0, 3) = -1.5*l1*l4*m4*sphi*dphi - 0.5*l1*l4*m4*sphi*dtheta1 - 1.5*l1*l4*m4*sphi*dtheta3;

    C(1, 0) = -1.5*l1*l2*(m2 + m3 + m4)*s12*dtheta1 - 0.5*l1*l2*(m2 + m3 + m4)*s12*dtheta2;
    C(1, 1) = -0.5*l1*l2*(m2 + m3 + m4)*s12*dtheta1 + 0.5*l1*l2*(m2 + m3 + m4)*s12*dtheta2;

    C(2, 0) = 0.5*l1*l4*m4*sphi*dphi + (-1.5*l1*l3*m3*s13 + 1.5*l1*l4*m4*sphi)*dtheta1
              + (-0.5*l1*l3*m3*s13 + 0.5*l1*l4*m4*sphi)*dtheta3;
    C(2, 2) = -0.5*l1*l4*m4*sphi*dphi + (-0.5*l1*l3*m3*s13 + 0.5*l1*l4*m4*sphi)*dtheta1
              + (0.5*l1*l3*m3*s13 - 0.5*l1*l4*m4*sphi)*dtheta3;
    C(2, 3) = -0.5*l1*l4*m4*sphi*dphi + 0.5*l1*l4*m4*sphi*dtheta1 - 0.5*l1*l4*m4*sphi*dtheta3;

    C(3, 0) = 0.5*l1*l4*m4*sphi*dphi + 1.5*l1*l4*m4*sphi*dtheta1 + 0.5*l1*l4*m4*sphi*dtheta3;
    C(3, 2) = -0.5*l1*l4*m4*sphi*dphi + 0.5*l1*l4*m4*sphi*dtheta1 - 0.5*l1*l4*m4*sphi*dtheta3;
    C(3, 3) = -0.5*l1*l4*m4*sphi*dphi + 0.5*l1*l4*m4*sphi*dtheta1 - 0.5*l1*l4*m4*sphi*dtheta3;

    // --- Gravity Vector G ---
    G(0) = g*l1*(m2 + m3 + m4)*cos(theta1) + g*m1*lambda1*cos(theta1);
    G(1) = g*l2*(m2 + m3 + m4)*cos(theta2) + g*m2*lambda2*cos(theta2);
    G(2) = g*l3*m3*cos(theta3) + g*l4*m4*cos(phi + theta3)
         + g*m3*lambda3*cos(theta3) + g*m4*lambda4*cos(phi + theta3);
    G(3) = g*l4*m4*cos(phi + theta3) + g*m4*lambda4*cos(phi + theta3);
}