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
        {50, 50, 50, 50});

    Eigen::Vector4d p(4);
    Eigen::Vector4d v(4);
    Eigen::Vector4d a(4);
    Eigen::Vector4d t(4);

    p.setOnes();
    v.setOnes();
    t.setZero();
    a.setZero();

    // std::cout << "A: " << std::endl;
    // std::cout << AD.discretizedA(p,v, 0.02) << std::endl;
    //
    //
    //


    // std::cout << "B: " << std::endl;
    // std::cout << AD.discretizedB(p,
    //     AD.A(p, v),
    //     AD.discretizedA(p, v, 0.02),
    //     0.02) << std::endl;
    // std::cout << AD.discretizedB(p, AD.A(p,v), AD.discretizedA(p, v, 0.02), 0.02) << std::endl;
    //
    // std::cout << "C: " << std::endl;
    // std::cout << AD.C() << std::endl;
    //
    // std::cout << "D: " << std::endl;
    // std::cout << AD.D() << std::endl;

    // MatrixExponentialReturnValue<>

    // std::cout << AD.G(p) << std::endl;
    // std::cout << AD.C(p, v) * v << std::endl;
    // std::cout << AD.M(p).inverse() << std::endl;
    // std::cout << AD.M(p).inverse() * (t - (AD.C(p,v)*v) - AD.G(p)) << std::endl;
    // std::cout << AD.getAccelerationsRPSS(t, p, v);

    int i = 0;
    while (i < 25) {
        a = AD.getAccelerationsRPSS(AD.G(p), p, v);
        std::cout << "A: " << a << std::endl;
        v = v + (0.2 * a);
        std::cout << "V: " << v << std::endl;
        p = p + (0.2 * v);
        std::cout << "P: " << p << std::endl;
        i++;
    }

    // std::cout << AD.M(p).exp() << std::endl;

    return 0;
}