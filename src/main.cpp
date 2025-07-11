#include <Eigen/Dense>
#include <iostream>
#include "Arm/ArmDynamics.h"

using namespace Eigen;

int main() {
    // Matrix2d mat;
    // mat << 1, 2,
    //        3, 4;
    // std::cout << mat << std::endl;

    ArmDynamics AD = ArmDynamics(
        {1, 1, 1, 1},
        {1, 1, 1, 1},
        {1, 1, 1, 1});

    Eigen::Vector4d p(4);
    Eigen::Vector4d v(4);
    Eigen::Vector4d a(4);
    Eigen::Vector4d t(4);

    p.setOnes();
    v.setZero();
    t.setZero();
    a.setZero();

    // std::cout << AD.M(p).inverse() * AD.M(p);

    while (true) {
        a == AD.getAccelerationsRPSS(t, p, v);
        std::cout << "A: " << a << std::endl;
        v = v + (0.01 * a);
        std::cout << "V: " << v << std::endl;
        p = p + (0.01 * v);
        std::cout << "P: " << p << std::endl;
    }

    return 0;
}