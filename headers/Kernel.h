#pragma once
#include <algorithm>

#include "Grid.h"

double CubicSplineKernel(const Eigen::Vector3d& r);
Eigen::Vector3d CubicSplineKernelGradient(const Eigen::Vector3d& r);

double CubicSplineKernel2D(Eigen::Vector2d r);
Eigen::Vector2d CubicSplineKernelGradient2D(Eigen::Vector2d r);

double CohesionSpline2D(Eigen::Vector2d r);
