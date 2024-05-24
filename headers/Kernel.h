#pragma once
#include <algorithm>

#include "Grid.h"

float CubicSplineKernel(const Eigen::Vector3f& r);
Eigen::Vector3f CubicSplineKernelGradient(const Eigen::Vector3f& r);

float CubicSplineKernel2D(Eigen::Vector2f r);
Eigen::Vector2f CubicSplineKernelGradient2D(Eigen::Vector2f r);

float CohesionSpline2D(Eigen::Vector2f r);
