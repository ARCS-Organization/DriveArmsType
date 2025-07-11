cmake_minimum_required(VERSION 3.28)
project(DriveArmsType)

set(CMAKE_CXX_STANDARD 17)

include_directories(include)

add_executable(DriveArmsType src/main.cpp
        src/Arm/ArmDynamics.cpp
        src/Arm/ArmDynamics.h)
