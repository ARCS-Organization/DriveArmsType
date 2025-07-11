#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;

int main() {
    Matrix2d mat;
    mat << 1, 2,
           3, 4;
    std::cout << mat << std::endl;
    return 0;
}