#pragma once
#include <algorithm>

#include "Grid.h"

float CubicSplineKernel(const Eigen::Vector3f& r);
Eigen::Vector3f CubicSplineKernelGradient(const Eigen::Vector3f& r);
